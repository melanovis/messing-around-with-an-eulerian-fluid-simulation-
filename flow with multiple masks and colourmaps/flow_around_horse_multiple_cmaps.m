format compact
clear
clc
%close all
clf reset

warning('off','all')

render_thing = false;
save_sim = true;

max_particle_quantity = 1e7;
micro_turbulance_factor_default = 1e-6;
grid_spacing = 1;
emission_turbulance_factor = 1; % 1 for complete intra pixel rand on the source mask
frame_rem_factor = 15; %higher is more simtime between frames
fps = 60;
shock_frames = [1]; %iters when shocks occur
maxiters = fps * frame_rem_factor * 30;
texteffect_scale_factor = 1;

png_name = "img_write.png";

collider_mask_raw = imread("encounter_mask.png");
collider_mask = single(squeeze(collider_mask_raw(:,:,1)));
collider_mask = flipud(collider_mask);
collider_mask(collider_mask<100)=0; %convert to logical

mask_extra_1_raw = imread("shoulder_mask.png");
mask_extra_1 = single(squeeze(mask_extra_1_raw(:,:,1)));
mask_extra_1 = flipud(mask_extra_1);
mask_extra_1 = logical(mask_extra_1);

mask_extra_2_raw = imread("horn_out_mask.png");
mask_extra_2 = single(squeeze(mask_extra_2_raw(:,:,1)));
mask_extra_2 = flipud(mask_extra_2);
mask_extra_2 = logical(mask_extra_2);

%alternate sourcemask
alternate_sourcemask = mask_extra_1;

%perforating collider mask
collider_mask_perf = collider_mask;
for n=1:height(collider_mask)
    for m=1:width(collider_mask)
        if collider_mask(n,m) && rand() > 0.01
            collider_mask_perf(n,m) = 1;
        end
    end
end

augmented_abyss = abyss;
augmented_abyss = augmented_abyss(round(linspace(1,height(augmented_abyss),5)),:);
cmap_1 = interp1( [0,1e-7,1e-6,1e-5,1e-3,1] , [[0, 0, 0]; augmented_abyss], linspace(0, 1, 1e3));
cmap_2 = interp1( [0,0.03,0.06,0.1,0.2,1] , [[0, 0, 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
%interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3)); 

cmap_series(:,:,1) = single(cmap_1);
cmap_series(:,:,2) = single(cmap_2);

text_scale_series = [
1
1
];

load("tiles_formatted.mat")

filename = "frame_store.mat";
if ~exist(filename, 'file')
    frame_colourised = [];
    save(filename,"frame_colourised")
end

filename_storesim = "storesim.mat";
if ~exist(filename_storesim, 'file')
    iter = 1;
    save(filename_storesim,"iter")
end

textcells_message = floor(rand(1,1e6)*16); 
textcells_message = string(dec2hex(textcells_message)).';

scene_scale = [height(collider_mask),width(collider_mask)]; %height, width

smoothing_filter=[
0.7, 0.7, 0.7 
0.7, 0.7, 0.7
0.7, 0.7, 0.7
];
smoothing_filter = smoothing_filter./numel(smoothing_filter);
smoothing_filter = smoothing_filter.* (1/sum(sum(smoothing_filter)));

[grid_x, grid_y] = meshgrid(1:scene_scale(2), 1:scene_scale(1));

[pressure_field,v_x,v_y] = deal( zeros(scene_scale(1), scene_scale(2)) );

%solid_mask = zeros(scene_scale(1), scene_scale(2));
solid_mask = collider_mask_perf > 125;
solid_mask = solid_mask + mask_extra_2;
solid_mask = logical(solid_mask);

if ~exist("solid_mask_store.mat", 'file') %this takes so long we need to be saving as we go
    save("solid_mask_store.mat","solid_mask")
else
    load("solid_mask_store.mat")
end

particle_inflow_mask = zeros(scene_scale(1), scene_scale(2));
particle_inflow_mask(3:end-2,end-2) = 1;
particle_inflow_mask = logical(particle_inflow_mask);
particle_inflow_mask = particle_inflow_mask + alternate_sourcemask;

outflow_bounds_x = [min(min(grid_x))+2,max(max(grid_x))-2];
outflow_bounds_y = [min(min(grid_y))+2,max(max(grid_y))-2];

%pressure masks
pressure_in_mask = zeros(scene_scale(1), scene_scale(2));
pressure_out_mask = zeros(scene_scale(1), scene_scale(2));
pressure_in_mask(3:end-2,end) = 1;
pressure_out_mask(end*1/3:end*2/3,1) = 1;
pressure_in_mask = logical(pressure_in_mask);
pressure_out_mask = logical(pressure_out_mask);
pressure_in_mask = flipud(pressure_in_mask);
pressure_out_mask = flipud(pressure_out_mask);

%regarding particles
particle_index = 1;
for n=1:height(particle_inflow_mask)
    for m=1:width(particle_inflow_mask)
        
        if particle_inflow_mask(n,m)
            particle_sourcemask_y(particle_index,1) = n;
            particle_sourcemask_x(particle_index,1) = m;
            
            if ~alternate_sourcemask(n,m)
                particle_source_ID(particle_index) = 1; %colourmap 1
            else
                particle_source_ID(particle_index) = 2; %colourmap 2
            end
                
            particle_index = particle_index+1;
        end
    end
end
particle_source_IDs = unique(particle_source_ID);

particlelist_x = [];
particlelist_y = [];
particle_ID_list = [];

scene_height = height(pressure_field);
scene_width = width(pressure_field);

if render_thing
    v = VideoWriter("fluid_horsemask_test", 'MPEG-4');
    v.FrameRate = fps;
    open(v);
end


iter = 1;
frame = 1;
while iter <= maxiters
    
    if exist(filename_storesim, 'file') && iter == 1
        load(filename_storesim);
    end

    v_x(solid_mask) = 0;
    v_y(solid_mask) = 0;
    
    v_divergence = divergence(v_x, v_y)/2;

    p_p = ones([size(v_x)]);
    d_p = ones([size(v_x)]);
    
    p_tmpfield = zeros(scene_scale(1), scene_scale(2));

    if iter > 5
        solve_iter_margin = 5e-3;
    else
        solve_iter_margin = 1e-3;
    end

    pressure_solve_iters = 1;
    while max(max(d_p)) > solve_iter_margin
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
        pressure_field(pressure_out_mask) = 0;
        pressure_field(pressure_in_mask) = 2;

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

    if ismember(iter,shock_frames)
        micro_turbulance_factor = 3;
    else
        micro_turbulance_factor = micro_turbulance_factor_default;
    end

    v_x = v_x + (rand(size(v_x))-0.5).*2*micro_turbulance_factor;
    v_y = v_y + (rand(size(v_y))-0.5).*2*micro_turbulance_factor;

    for n=1:6 %how many particles are we adding?
        particlelist_x = [particle_sourcemask_x + rand(height(particle_sourcemask_x),1)*emission_turbulance_factor; particlelist_x];
        particlelist_y = [particle_sourcemask_y + rand(height(particle_sourcemask_y),1)*emission_turbulance_factor; particlelist_y];
        particle_ID_list = [particle_source_ID, particle_ID_list];
    end

    [particlelist_x, particlelist_y] = RK4_step(particlelist_x, particlelist_y, v_x, v_y, 1); %forward advection

    %cull particles on boundary
    particle_boundary_mask = [particlelist_x < outflow_bounds_x(1)] + [particlelist_x > outflow_bounds_x(2)] + [particlelist_y < outflow_bounds_y(1)] + [particlelist_y > outflow_bounds_y(2)];
    particle_boundary_mask = logical(particle_boundary_mask);
    particlelist_x(particle_boundary_mask) = [];
    particlelist_y(particle_boundary_mask) = [];
    particle_ID_list(particle_boundary_mask) = [];

    %particles intersecting with solid
    ind_colliding = find_colliding_particles(particlelist_x, particlelist_y, solid_mask);
    particlelist_x(ind_colliding) = [];
    particlelist_y(ind_colliding) = [];
    particle_ID_list(ind_colliding) = [];

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
    particle_ID_list(particle_ind_remove) = [];

    
    if rem(iter,frame_rem_factor) == 0

        frame_colour_total = uint8(zeros(size(pressure_field)));

        for ind_ID = 1:length(particle_source_IDs)

            texteffect_scale_factor = text_scale_series(ind_ID);

            particlelist_x_spec = particlelist_x(particle_ID_list == ind_ID);
            particlelist_y_spec = particlelist_y(particle_ID_list == ind_ID);

            %build density map
            particle_densitymap = particlelist_to_densitymap(particlelist_x_spec, particlelist_y_spec, [size(pressure_field)]);
            particle_densitymap = particle_densitymap./max(max(particle_densitymap));
    
            %smooth out density map
            smooth_densitymap = smooth_map(particle_densitymap,smoothing_filter);
            
            plot_field = smooth_densitymap;
        
            plot_field = plot_field./max(max(plot_field));
    
            plot_solidmask = single(solid_mask);

            textgrid_size = [round(scene_height/texteffect_scale_factor),round(scene_width/texteffect_scale_factor)];
            plotfield_textcells = imresize(plot_field,[textgrid_size(1),textgrid_size(2)]);
            plotfield_textcells(plotfield_textcells<1e-3) = 0;
            
            solidmask_scaledown = imresize(plot_solidmask,[textgrid_size(1),textgrid_size(2)]);
            solidmask_scaledown(solidmask_scaledown<1e-3) = 0;
    
            cmap = squeeze(cmap_series(:,:,ind_ID));

            frame_colourised = [];
            if texteffect_scale_factor <= 1 
                frame_colourised = build_coloured_videoframe(plot_field, solid_mask, cmap, false);
                return_tilemap = plot_field; %used for overwriting prev particlemap if necessary
            else
                textcell_order = build_order_map_with_solidmask(plotfield_textcells,solidmask_scaledown);
                size_original = size(plot_field);
                [textcelled_flow, return_tilemap] = implement_chars_on_map(textcell_order,textcells_message,textgrid_size,texteffect_scale_factor,size_original,plotfield_textcells,textcell_order,char_tiles,tile_names);
                textcelled_flow = textcelled_flow./max(max(textcelled_flow));
                frame_colourised = build_coloured_videoframe(textcelled_flow, solid_mask, cmap, false);
            end

            %bad way of masking out prev flows but oh well
            if ind_ID>1
                for n=1:height(frame_colourised)
                    for m=1:width(frame_colourised)
                        if sum(frame_colour_total(n,m,:)) > 0 && return_tilemap(n,m) > 0
                            frame_colour_total(n,m,:) = 0;
                        end
                    end
                end
            end
            frame_colour_total = frame_colour_total + uint8(frame_colourised);

        end

        scatter(nan,nan)
        hold on
        imagesc(frame_colour_total)
        imagesc(solid_mask, alphadata = solid_mask)
        axis equal tight
        set(gca,"ColorScale","log")
        hold off

        set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

        drawnow()
        iter


        frame_series(:,:,frame) = single(plot_field);
        frame = frame+1;
        if save_sim
            save(filename,"frame_colourised","frame_colour_total","solid_mask","plot_field");
            imwrite(flipud(frame_colour_total),png_name)
        end

        if save_sim
            save(filename_storesim,"pressure_field","particlelist_x","particlelist_y","particle_ID_list","iter","v_x","v_y")
        end
    end

    iter = iter+1;

end

if render_thing
    close(v);
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);



function return_frame = greyscale_img_solidmask(solid_mask)
    return_frame = uint8(zeros([size(solid_mask),3]));
    for n=1:3
        return_frame(:,:,n) = solid_mask.*255;
    end
end

function return_frame = build_coloured_videoframe(plot_field, solid_mask, cmap, add_solidmask)
    scene_scale = size(plot_field);
    plot_field_norm = plot_field;
    scaling_map = flip( round(height(cmap) - logspace(log10(1),log10(height(cmap)), height(cmap) )) + 1 );
    %scaling_map = 1:height(cmap);
    plot_field_norm = ceil( plot_field_norm.*height(cmap) );
    plot_field_norm(plot_field_norm==0) = 1;
    for n=1:3
        cmap_channel_spec = cmap(scaling_map(plot_field_norm),n) .* 255;
        %cmap_channel_spec = cmap(plot_field_norm,n) .* 255;
        field_colour_spec = reshape(cmap_channel_spec,[scene_scale]);
        if add_solidmask
            field_colour_spec(solid_mask) = 255;
        end
        return_frame(:,:,n) = single(field_colour_spec);
    end
end


function [return_map,return_tilemap] = implement_chars_on_map(order_map,message,textgrid_size,scale_factor,size_original,plotfield_textcells,textcell_order,char_tiles,tile_names)
    
    return_map = zeros(size_original);
    return_tilemap = zeros(size_original);
    order_linear = reshape(order_map(order_map~=0),1,[]);

    if length(message) < length(order_linear)
        message = [message, repelem("0",length(order_linear) - length(message))];
    end

    ind_order = 1;

    for y=height(order_map):-1:1
        for x=1:width(order_map)
            if order_map(y,x) ~= 0
                %fprintf("------\n")
                re_x1 = (x-1)*scale_factor+1;
                re_y1 = (y-1)*scale_factor+1;
                re_x2 = re_x1+scale_factor-1;
                re_y2 = re_y1+scale_factor-1;
                
                %return_map(re_y1:re_y2,re_x1:re_x2) = order_linear(ind_order);

                %what character are we reading?
                character_read = message(ind_order);
                char_tile_spec = getfield(char_tiles,"tile_"+character_read);

                char_tile_spec = imresize(char_tile_spec,[re_y2-re_y1+1,nan]);
                char_tile_spec = flipud(abs(char_tile_spec));
                char_tile_spec(char_tile_spec<0.1) = 0;

                return_tilemap(re_y1:re_y2,re_x1:re_x2) = 1;

                if rand() < 1e-4
                    char_tile_spec = circshift(char_tile_spec,[0,1]);
                    char_tile_spec(:,1) = 0;
                end
                return_map(re_y1:re_y2,re_x1:re_x1+width(char_tile_spec)-1) = char_tile_spec .* plotfield_textcells(y,x);

                ind_order = ind_order+1;
            end
        end
    end
    return_map = return_map(1:size_original(1),1:size_original(2));
end

function order_map = build_order_map_with_solidmask(grid,solid_mask)
    ind=1;
    order_map = zeros(size(grid));
    for m=1:width(grid)
        for n=1:height(grid)
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

    %smoothing out single holes
    % hole_mask = find_single_holes(region_map);
    % fluidmap = ~region_map | hole_mask;
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