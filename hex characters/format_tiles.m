format compact
clear
clc
%close all
clf reset

namelist = ["0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"];

char_tiles = [];
tile_names = namelist;
tile_height = 200; %known

for n=1:length(namelist)
    char_raw = single( imread([namelist(n)+".png"]) );
    char_total = single( zeros(size(char_raw,1),size(char_raw,2)) );
    for m=1:size(char_raw,3)
        char_total = char_total + squeeze(char_raw(:,:,m));
    end
    char_total = 1 - char_total./max(max(char_total));
    char_tiles = setfield(char_tiles,"tile_"+namelist(n),char_total);
end

save("tiles_formatted.mat","tile_names","char_tiles","tile_height")