function [ ex,img,targetregion,imgCompleteFinal ] = extern_imextrapolate( img_filename,mode,load_img_postfix, resize_factor,im2bin_mode,erosion )
tic;
perc_for_n_pix = 0.067;
if ~exist('load_img_postfix','var')
    load_img_postfix = '';
end
if ~exist('resize_factor','var')
    resize_factor = 1;
end
[~,name,ext]=fileparts(img_filename);
img=imresize(imread([name,load_img_postfix,ext]),resize_factor);
img=s_canvasSize(img,size(img)+30);
targetregion = ~s_scale_mask(s_im2bin(img,rgb2lab(img),im2bin_mode),-erosion); %mask is the extrapolation target region
% img(s_ind2to3(targetregion))=255;
[r,c]=find(~targetregion);
r_l=max(r)-min(r);
c_l=max(c)-min(c);
n_px = round(mean([r_l*perc_for_n_pix,c_l*perc_for_n_pix]));
if strcmp(mode,'m')
    imgCompleteFinal = im_complete(img, targetregion);
elseif strcmp(mode,'c')
    imgCompleteFinal = sc_complete(img, targetregion);
else
    error('Unknown mode');
end
ex = imgCompleteFinal;
ex_bg_ind = s_ind2to3(~s_scale_mask(~targetregion,n_px));
ex(ex_bg_ind) = repmat(get_background_color(img),size(ex_bg_ind,1),1);
fprintf('Extrapolate on image %s took: %0.3f sec\n',img_filename,toc);
end

