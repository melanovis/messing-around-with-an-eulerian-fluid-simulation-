format compact
clear
clc
%close all
clf reset

filename = "img_write.png";

load("frame_store.mat")

hold on
grid on
axis tight equal
imagesc(frame_colour_total)
imwrite(frame_colour_total,filename)