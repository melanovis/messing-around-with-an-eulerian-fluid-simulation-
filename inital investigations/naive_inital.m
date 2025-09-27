format compact
clear
clc
%close all
clf reset

s = 500;
ar = 1.5;
grid_spacing = 1;

jacobi_mask = [
0 1 0 
1 0 1
0 1 0
]./4;

[grid_x, grid_y] = meshgrid(1:s*ar, 1:s);

pressure_field = zeros(s, s*ar);
v_x = zeros(s, s*ar);
v_y = zeros(s, s*ar);

solid_mask = zeros(s, s*ar); 
for n=10:20 %dummy solid mask
    solid_mask(50:100,n) = rand(1,abs(50-100)+1) < 0.4;
    %solid_mask(40:60,n) = rem(40:60,3);
end
solid_mask = logical(solid_mask);

[px, py] = meshgrid(10:10, 1:s);

iter = 0;

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
colormap(cmap)

while true

    tic

    v_x(solid_mask) = 0;
    v_y(solid_mask) = 0;

    v_divergence = divergence(v_x, v_y)/2;

    t1 = toc;

    for n = 1:10
        pressure_field = (conv2(pressure_field, jacobi_mask, 'same') - v_divergence);

        %boundary conditions
        pressure_field(1,1:end) = pressure_field(2,1:end);
        pressure_field(end,1:end) = pressure_field(end-1,1:end);
        
        pressure_field(1:end,end) = 0;

        pressure_field(3:end-2,1) = 1;
    end

    t2 = toc - t1;

    %compute velocity gradient and update velocities
    [dx, dy] = gradient(pressure_field);
    v_x(3:end-2, 3:end-2) = v_x(3:end-2, 3:end-2) - dx(3:end-2, 3:end-2);
    v_y(3:end-2, 3:end-2) = v_y(3:end-2, 3:end-2) - dy(3:end-2, 3:end-2);
    
    % dx = 0.5 * (pressure_field(3:end-2, 4:end-1) - pressure_field(3:end-2, 2:end-3)) ./ grid_spacing;
    % dy = 0.5 * (pressure_field(4:end-1, 3:end-2) - pressure_field(2:end-3, 3:end-2)) ./ grid_spacing;
    % v_x(3:end-2, 3:end-2) = v_x(3:end-2, 3:end-2) - dx;
    % v_y(3:end-2, 3:end-2) = v_y(3:end-2, 3:end-2) - dy;

    t3 = toc - t2;

    %advect
    [pv_x, pv_y] = RK4(grid_x, grid_y, v_x, v_y, -1);
    v_x = interp2(v_x, pv_x, pv_y, 'linear', 0);
    v_y = interp2(v_y, pv_x, pv_y, 'linear', 0);

    t4 = toc - t3;

    if rem(iter,20) == 0

        curl_field = curl(v_x, v_y);
        plot_field = abs(curl_field);

        plot_field(solid_mask) = 0;
        plot_field = plot_field(4:end-3,4:end-3);
        plot_field = plot_field./max(max(plot_field));
        
        surf(plot_field,EdgeColor="none");
        view([0,90])
        axis equal tight 
        %set(gca,"ColorScale","log")
        
        drawnow()
        iter
    end

    iter = iter+1;
end



function [x_new, y_new] = RK4(p_x, p_y, v_x, v_y, h)
   k1x = interp2(v_x, p_x, p_y, 'linear', 0);
   k1y = interp2(v_y, p_x, p_y, 'linear', 0);
   k2x = interp2(v_x, p_x + h/2 * k1x, p_y + h/2 * k1y, 'linear', 0);
   k2y = interp2(v_y, p_x + h/2 * k1x, p_y + h/2 * k1y, 'linear', 0);
   k3x = interp2(v_x, p_x + h/2 * k2x, p_y + h/2 * k2y, 'linear', 0);
   k3y = interp2(v_y, p_x + h/2 * k2x, p_y + h/2 * k2y, 'linear', 0);
   k4x = interp2(v_x, p_x + h * k3x, p_y + h * k3y, 'linear', 0);
   k4y = interp2(v_y, p_x + h * k3x, p_y + h * k3y, 'linear', 0);
   x_new = p_x + h/6 * (k1x + 2*k2x + 2*k3x + k4x);
   y_new = p_y + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
end