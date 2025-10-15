format compact
clear
clc
%close all
clf reset

load("storesim.mat");
load("solid_mask_store.mat")

cmap = interp1( [0,1.5e-2,2e-2,6e-2,0.15,1] , [[0, 0, 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3)); %not regular cmap

texteffect_scale_factor = 1;

smoothing_filter=[
0.7, 0.7, 0.7 
0.7, 0.7, 0.7
0.7, 0.7, 0.7
];

scene_height = height(pressure_field);
scene_width = width(pressure_field);

particle_densitymap = particlelist_to_densitymap(particlelist_x, particlelist_y, [size(pressure_field)]);
particle_densitymap = particle_densitymap./max(max(particle_densitymap));

%smooth out density map
smooth_densitymap = single( smooth_map(particle_densitymap,smoothing_filter) );

plot_field = smooth_densitymap;

plot_field = plot_field./max(max(plot_field));

plot_solidmask = single(solid_mask);

%building text map
textgrid_size = [round(scene_height/texteffect_scale_factor),round(scene_width/texteffect_scale_factor)];
plotfield_textcells = imresize(plot_field,[textgrid_size(1),textgrid_size(2)]);
plotfield_textcells(plotfield_textcells<1e-3) = 0;

solidmask_scaledown = imresize(plot_solidmask,[textgrid_size(1),textgrid_size(2)]);
solidmask_scaledown(solidmask_scaledown<1e-3) = 0;

if texteffect_scale_factor <= 1 
    frame_colourised = build_coloured_videoframe(plot_field, solid_mask, cmap, false);
else
    textcell_order = build_order_map_with_solidmask(plotfield_textcells,solidmask_scaledown);
    size_original = size(plot_field);
    textcelled_flow = implement_chars_on_map(textcell_order,textcells_message,textgrid_size,texteffect_scale_factor,size_original,plotfield_textcells,textcell_order,char_tiles,tile_names);
    textcelled_flow = textcelled_flow./max(max(textcelled_flow));
    frame_colourised = build_coloured_videoframe(textcelled_flow, solid_mask, cmap, false);
end

colormap(cmap)

scatter(nan,nan)
hold on
imagesc(frame_colourised)
imagesc(solid_mask, alphadata = solid_mask)
axis equal tight
set(gca,"ColorScale","log")
hold off

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

drawnow()



function return_frame = build_coloured_videoframe(plot_field, solid_mask, cmap, add_solidmask)
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
        if add_solidmask
            field_colour_spec(solid_mask) = 255;
        end
        return_frame(:,:,n) = uint8(field_colour_spec);
    end
end

function return_map = implement_chars_on_map(order_map,message,textgrid_size,scale_factor,size_original,plotfield_textcells,textcell_order,char_tiles,tile_names)
    
    return_map = zeros(size_original);
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

                % char_tile_spec = imresize(char_tile_spec,[re_y2-re_y1+1, re_x2-re_x1+1]); %make nan for no stretching
                % char_tile_spec = flipud(abs(char_tile_spec));
                % return_map(re_y1:re_y2,re_x1:re_x2) = char_tile_spec;

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

% function order_map = build_order_map(grid)
%     ind=1;
%     order_map = zeros(size(grid));
%     for n=1:height(grid)
%         for m=width(grid):-1:1
%             if grid(n,m)>0
%                 order_map(n,m) = ind;
%                 ind=ind+1;
%             end
%         end
%     end
% end

function smooth_densitymap = smooth_map(particle_densitymap,smoothing_filter)
    region_map = logical(particle_densitymap == 0);
    smooth_densitymap = conv2(particle_densitymap, smoothing_filter, 'same');

    fluidmap = ~region_map;

    smooth_densitymap(~fluidmap) = 0;
end

function return_densitymap = particlelist_to_densitymap(particlelist_x, particlelist_y,scene_scale)
    particlelist_x = round(particlelist_x);
    particlelist_y = round(particlelist_y);
    ind_linear = sub2ind([scene_scale], particlelist_y, particlelist_x);
    return_densitymap = accumarray(ind_linear, 1, [scene_scale(1) * scene_scale(2), 1]);
    return_densitymap = reshape(return_densitymap, [scene_scale]);
end