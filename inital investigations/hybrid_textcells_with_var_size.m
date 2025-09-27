format compact
clear
clc
%close all
clf reset


render_test = true;

scene_scale = [1e3,1e3]; %height, width
max_particle_quantity = 1e6;
micro_turbulance_factor_default = 1e-7;
grid_spacing = 1;
emission_turbulance_factor = 1; % 1 for complete intra pixel rand on the source mask
frame_rem_factor = 10; %higher is more simtime between frames
fps = 60;

load("tiles_formatted.mat")

textcells_message = floor(rand(1,1e6)*16); %rand for now, real later
textcells_message = string(dec2hex(textcells_message)).';

maxiters = fps * frame_rem_factor * 30;

scale_series_unique = flip( unique([20:-2:1,1]) );
testing_scale_series = ones([1,maxiters]);
testing_scale_series(1:end/2) = round(imresize(scale_series_unique,[1,maxiters/2],"nearest")); %evenly space out scaled character scale series

smoothing_filter=[
0.7, 0.7, 0.7 
0.7, 0.7, 0.7
0.7, 0.7, 0.7
];
smoothing_filter = smoothing_filter./numel(smoothing_filter);
smoothing_filter = smoothing_filter.* (1/sum(sum(smoothing_filter)));

[grid_x, grid_y] = meshgrid(1:scene_scale(2), 1:scene_scale(1));

[pressure_field,v_x,v_y] = deal( zeros(scene_scale(1), scene_scale(2)) );

solid_mask = zeros(scene_scale(1), scene_scale(2));
for n=50:70
    solid_mask(100:900,n) = rand(1,abs(100-900)+1) < 0.05;
end
solid_mask = logical(solid_mask);

particle_inflow_mask = zeros(scene_scale(1), scene_scale(2));
particle_inflow_mask(150:850,3) = 1;
particle_inflow_mask = logical(particle_inflow_mask);

outflow_bounds_x = [min(min(grid_x))+2,max(max(grid_x))-2];
outflow_bounds_y = [min(min(grid_y))+2,max(max(grid_y))-2];

%regarding particles
particle_index = 1;
for n=1:height(particle_inflow_mask)
    for m=1:width(particle_inflow_mask)
        if particle_inflow_mask(n,m)
            particle_sourcemask_y(particle_index,1) = n;
            particle_sourcemask_x(particle_index,1) = m;
            particle_index = particle_index+1;
        end
    end
end

iter = 1;

particlelist_x = [];
particlelist_y = [];

if render_test
    v = VideoWriter("fluid_textcell_test", 'MPEG-4');
    v.FrameRate = fps;
    open(v);
end

scene_height = height(pressure_field);
scene_width = width(pressure_field);

while iter <= maxiters

    texteffect_scale_factor = testing_scale_series(iter);

    v_x(solid_mask) = 0;
    v_y(solid_mask) = 0;
    
    v_divergence = divergence(v_x, v_y)/2;

    p_p = ones([size(v_x)]);
    d_p = ones([size(v_x)]);
    
    p_tmpfield = zeros(scene_scale(1), scene_scale(2));

    pressure_solve_iters = 1;
    while max(max(d_p)) > 1e-3

        p_tmpfield_prev = pressure_field;

        for n = 2:scene_height-1
            for m = 2:scene_width-1
                p_tmpfield(n,m) = 0.25 * (p_tmpfield_prev(n+1,m) + p_tmpfield_prev(n-1,m) + ...
                p_tmpfield_prev(n,m+1) + p_tmpfield_prev(n,m-1) - ...
                v_divergence(n,m));
            end
        end
        pressure_field = p_tmpfield;

        %boundary conditions
        pressure_field(1,1:end) = pressure_field(2,1:end);
        pressure_field(end,1:end) = pressure_field(end-1,1:end);

        %flow constraints
        pressure_field(1:end/2,end) = 0;
        pressure_field(3:end-2,1) = 1;

        d_p = abs(pressure_field - p_p);
        p_p = pressure_field;
        pressure_solve_iters = pressure_solve_iters+1;
    end
    
    dx = 0.5 * (pressure_field(3:end-2, 4:end-1) - pressure_field(3:end-2, 2:end-3)) ./ grid_spacing;
    dy = 0.5 * (pressure_field(4:end-1, 3:end-2) - pressure_field(2:end-3, 3:end-2)) ./ grid_spacing;
    v_x(3:end-2, 3:end-2) = v_x(3:end-2, 3:end-2) - dx;
    v_y(3:end-2, 3:end-2) = v_y(3:end-2, 3:end-2) - dy;

    [pv_x, pv_y] = RK4_step(grid_x, grid_y, v_x, v_y, -1); %backward advection
    v_x = interp2(v_x, pv_x, pv_y, 'linear', 0);
    v_y = interp2(v_y, pv_x, pv_y, 'linear', 0);

    %micro_turbulance_factor = micro_turbulance_factor_default;
    if rem(iter,3000) == 0
        micro_turbulance_factor = 3; %sudden shock
    else
        micro_turbulance_factor = micro_turbulance_factor_default;
    end

    v_x = v_x + (rand(size(v_x))-0.5).*2*micro_turbulance_factor;
    v_y = v_y + (rand(size(v_y))-0.5).*2*micro_turbulance_factor;

    for n=1:2 %how many particles are we adding?
        particlelist_x = [particle_sourcemask_x + rand(height(particle_sourcemask_x),1)*emission_turbulance_factor; particlelist_x];
        particlelist_y = [particle_sourcemask_y + rand(height(particle_sourcemask_y),1)*emission_turbulance_factor; particlelist_y];
    end

    [particlelist_x, particlelist_y] = RK4_step(particlelist_x, particlelist_y, v_x, v_y, 1); %forward advection

    %cull particles on boundary
    particle_boundary_mask = [particlelist_x < outflow_bounds_x(1)] + [particlelist_x > outflow_bounds_x(2)] + [particlelist_y < outflow_bounds_y(1)] + [particlelist_y > outflow_bounds_y(2)];
    particle_boundary_mask = logical(particle_boundary_mask);
    particlelist_x(particle_boundary_mask) = [];
    particlelist_y(particle_boundary_mask) = [];

    %particles intersecting with solid
    ind_colliding = find_colliding_particles(particlelist_x, particlelist_y, solid_mask);
    particlelist_x(ind_colliding) = [];
    particlelist_y(ind_colliding) = [];

    %remove when too many particles
    particle_ind_remove = [];
    if length(particlelist_x) > max_particle_quantity
        cull_quantity = length(particlelist_x) - max_particle_quantity;
        while length(particle_ind_remove) < cull_quantity
            particle_ind_remove = [particle_ind_remove, ceil( rand(1,length(cull_quantity))*length(particlelist_x) )];
            particle_ind_remove = unique(particle_ind_remove);
        end
    end
    particlelist_x(particle_ind_remove) = [];
    particlelist_y(particle_ind_remove) = [];

    
    if rem(iter,frame_rem_factor) == 0

        %build density map
        particle_densitymap = particlelist_to_densitymap(particlelist_x, particlelist_y, [size(pressure_field)]);
        particle_densitymap = particle_densitymap./max(max(particle_densitymap));

        %smooth out density map
        smooth_densitymap = smooth_map(particle_densitymap,smoothing_filter);

        plot_field = smooth_densitymap;
    
        plot_field = plot_field./max(max(plot_field));

        plot_solidmask = single(solid_mask);

        %building text map
        textgrid_size = [round(scene_height/texteffect_scale_factor),round(scene_width/texteffect_scale_factor)];
        plotfield_textcells = imresize(plot_field,[textgrid_size(1),textgrid_size(2)]);
        plotfield_textcells(plotfield_textcells<1e-3) = 0;
        
        solidmask_scaledown = imresize(plot_solidmask,[textgrid_size(1),textgrid_size(2)]);
        solidmask_scaledown(solidmask_scaledown<1e-3) = 0;

        textcell_order = build_order_map_with_solidmask(plotfield_textcells,solidmask_scaledown);
        %textcell_order = build_order_map(plotfield_textcells);
        size_original = size(plot_field);
        textcelled_flow = implement_chars_on_map(textcell_order,textcells_message,textgrid_size,texteffect_scale_factor,size_original,plotfield_textcells,textcell_order,char_tiles,tile_names);
        textcelled_flow = textcelled_flow./max(max(textcelled_flow));
                    
        cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
        colormap(cmap)
        
        if texteffect_scale_factor <= 1 
            frame_colourised = build_coloured_videoframe(plot_field, solid_mask, cmap);
        else
            frame_colourised = build_coloured_videoframe(textcelled_flow, solid_mask, cmap);
        end

        scatter(nan,nan)
        hold on
        imagesc(frame_colourised)
        axis equal tight
        set(gca,"ColorScale","log")
        hold off

        set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

        drawnow()
        iter

        if render_test
            writeVideo(v,flipud(frame_colourised))
        end
    end

    iter = iter+1;
end

if render_test
    close(v);
end


function return_frame = build_coloured_videoframe(plot_field, solid_mask, cmap)
    scene_scale = size(plot_field);
    plot_field_norm = plot_field;
    %colourising
    logscale_map = flip( round(height(cmap) - logspace(log10(1),log10(height(cmap)), height(cmap) )) + 1 );
    plot_field_norm = ceil( plot_field_norm.*height(cmap) );
    plot_field_norm(plot_field_norm==0) = 1;    
    for n=1:3
        cmap_channel_spec = cmap(logscale_map(plot_field_norm),n) .* 255;
        %cmap_channel_spec = cmap(plot_field_norm,n) .* 255;
        field_colour_spec = reshape(cmap_channel_spec,[scene_scale]);
        field_colour_spec(solid_mask) = 255;
        return_frame(:,:,n) = uint8(field_colour_spec);
    end
end

function return_map = implement_chars_on_map(order_map,message,textgrid_size,scale_factor,size_original,plotfield_textcells,textcell_order,char_tiles,tile_names)
    
    return_map = zeros(size_original);
    order_linear = reshape(order_map(order_map~=0),1,[]);

    if length(message) < length(order_linear)
        error("message has to be longer than grid order, fix me later")
    end

    ind_order = 1;

    for x=1:width(order_map)
        for y=1:height(order_map)
            if order_map(y,x) ~= 0
                %fprintf("------\n")
                re_x1 = (x-1)*scale_factor+1;
                re_y1 = (y-1)*scale_factor+1;
                re_x2 = re_x1+scale_factor-1;
                re_y2 = re_y1+scale_factor-1;

                %what character are we reading?
                character_read = message(ind_order);
                char_tile_spec = getfield(char_tiles,"tile_"+character_read);
                char_tile_spec = imresize(char_tile_spec,[re_y2-re_y1+1, nan]); %make nan for no stretching
                char_tile_spec = flipud(abs(char_tile_spec));
                char_tile_spec(char_tile_spec<0.1) = 0;
                if rand() < 1e-4 && character_read~="D" %some random jitter (character D may clip)
                    char_tile_spec = circshift(char_tile_spec,[0,1]);
                    char_tile_spec(:,1) = 0;
                end
                return_map(re_y1:re_y2,re_x1:re_x1+width(char_tile_spec)-1) = char_tile_spec .* plotfield_textcells(y,x);

                ind_order = ind_order+1;
            end
        end
    end
    return_map = return_map(1:size_original(2),1:size_original(1));

end

function order_map = build_order_map_with_solidmask(grid,solid_mask) %solid map is excluded
    ind=1;
    order_map = zeros(size(grid));
    for n=1:height(grid)
        for m=width(grid):-1:1
            if grid(n,m) > 0 && solid_mask(n,m) == 0
                order_map(n,m) = ind;
                ind=ind+1;
            end
        end
    end
end

function order_map = build_order_map(grid)
    ind=1;
    order_map = zeros(size(grid));
    for n=1:height(grid)
        for m=width(grid):-1:1
            if grid(n,m)>0
                order_map(n,m) = ind;
                ind=ind+1;
            end
        end
    end
end

function smooth_densitymap = smooth_map(particle_densitymap,smoothing_filter)
    region_map = logical(particle_densitymap == 0);
    smooth_densitymap = conv2(particle_densitymap, smoothing_filter, 'same');
    fluidmap = ~region_map;
    smooth_densitymap(~fluidmap) = 0;
end

function hole_mask = find_single_holes(mask)
    [mask_up, mask_down, mask_left, mask_right] = deal(false(size(mask)));
    mask_up(1:end-1, :) = ~mask(2:end, :);
    mask_down(2:end, :) = ~mask(1:end-1, :);
    mask_left(:, 1:end-1) = ~mask(:, 2:end);
    mask_right(:, 2:end) = ~mask(:, 1:end-1);
    hole_mask = mask & mask_up & mask_down & mask_left & mask_right;
end

function ind_colliding = find_colliding_particles(particlelist_x, particlelist_y, solid_mask)
    particle_x_closest = round(particlelist_x);
    particle_y_closest = round(particlelist_y);
    [solid_coord_y,solid_coord_x] = find(solid_mask); 
    [is_colliding, ~] = ismember([particle_x_closest,particle_y_closest], [solid_coord_x,solid_coord_y], 'rows');
    ind_colliding = find(is_colliding);
end

function return_densitymap = particlelist_to_densitymap(particlelist_x, particlelist_y,scene_scale)
    particlelist_x = round(particlelist_x);
    particlelist_y = round(particlelist_y);
    ind_linear = sub2ind([scene_scale], particlelist_y, particlelist_x);
    return_densitymap = accumarray(ind_linear, 1, [scene_scale(1) * scene_scale(2), 1]);
    return_densitymap = reshape(return_densitymap, [scene_scale]);
end

function [x_new, y_new] = RK4_step(p_x, p_y, v_x, v_y, h)
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