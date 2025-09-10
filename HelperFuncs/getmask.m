function [ mask ] = getmask( pics_dir, filename )
if exist(pics_dir,'dir')
    addpath(pics_dir);
else
    error('pics directory not found');
end
im2bin_mode='other'; 
addpath(genpath(pwd));
m=imread(filename);
mask = s_im2bin(m,rgb2lab(m),im2bin_mode);
end

