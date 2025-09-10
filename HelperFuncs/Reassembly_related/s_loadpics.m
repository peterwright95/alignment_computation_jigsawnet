function [ input,fr_info ] = s_loadpics( pics, options)
EXTRAPOLATED_POSTFIX = '_ext';
if ~exist('options','var')
    options = add_missing_options([]);
else
    options = add_missing_options(options);
end

SCALE_MASK = -1;
INIT_PAD_SIZE = 10;
im2bin_mode='other'; %other/palette

n_pics = numel(pics);
input = cell(n_pics,1);
fr_info = cell(n_pics,1);
for i=1:n_pics
    filename = pics{i};
    [~,name,ext] = fileparts(filename);
    filename = [name ext];
    filenameext = [name EXTRAPOLATED_POSTFIX ext];
    
    if options.is_use_rand_frag_rot
        if isfield(options, 'rand_angles')
            rand_angle = options.rand_angles(i);
        else
            rand_angle = round(rand*359);
        end        
    else
        rand_angle = 0;
    end
    input{i}.rand_angle = rand_angle;
    im = imread(filename);
    im = s_canvasSize(im,size(im)+INIT_PAD_SIZE);
    bgc = get_background_color(im);
    after_rot_partial_bg_i = s_ind2to3(~imrotate(true(size(im,1),size(im,2)),rand_angle));
    im = imrotate(im,rand_angle);
    im(after_rot_partial_bg_i) = repmat(bgc,size(after_rot_partial_bg_i,1),1);
    input{i}.rgb = im;
    [input{i}.origsz(1),input{i}.origsz(2),~] = size(input{i}.rgb);
    
    im_ex = imread(filenameext);
    im_ex = s_canvasSize(im_ex,size(im_ex)+INIT_PAD_SIZE);
    ex_bgc = get_background_color(im);
    after_rot_partial_bg_i = s_ind2to3(~imrotate(true(size(im_ex,1),size(im_ex,2)),rand_angle));
    im_ex = imrotate(im_ex,rand_angle);
    im_ex(after_rot_partial_bg_i) = repmat(ex_bgc,size(after_rot_partial_bg_i,1),1);
    input{i}.rgb_ex = im_ex;
    [input{i}.origsz_ex(1),input{i}.origsz_ex(2),~] = size(input{i}.rgb_ex);
    
    assert(isequal(input{i}.origsz,input{i}.origsz_ex),'regular and extrapolated fragment images should be of equal size');
    
    input{i}.lab = rgb2lab(input{i}.rgb);
    input{i}.lab_ex = rgb2lab(input{i}.rgb_ex);
    input{i}.mask = s_im2bin(input{i}.rgb,input{i}.lab,im2bin_mode);
    input{i}.mask = s_scale_mask(input{i}.mask,SCALE_MASK);
    [input{i}.mask_Si,input{i}.mask_Mi] = s_reorder_mask(input{i}.mask);%input{i}.mask_Si
    input{i}.mask_ex = s_im2bin(input{i}.rgb_ex,input{i}.lab_ex,im2bin_mode);
    input{i}.mask_ex = s_scale_mask(input{i}.mask_ex,SCALE_MASK);
    input{i}.poly = s_mask2poly(input{i}.mask,options.reduce_poly_length_to);
    input{i}.poly_ex = s_mask2poly(input{i}.mask_ex,options.reduce_poly_length_to);
    palette_sz = [3,4,6,8,10];
    [input{i}.grad_ex,input{i}.gx_ex,input{i}.gy_ex] = s_palette_grad(input{i}.rgb_ex,input{i}.lab_ex,palette_sz,'sobel',options.is_use_dist_colors_4grad,options.is_smooth_before_grad);
    if options.is_grad_dir_on_plt_output
        [input{i}.gx_ex,input{i}.gy_ex]=imgradientxy(input{i}.grad_ex);
    end
    input{i}.gx = input{i}.gx_ex;
    input{i}.gy = input{i}.gy_ex;
    input{i}.grad = input{i}.grad_ex;
    
    input{i}.gx(~input{i}.mask) = 0;
    input{i}.gy(~input{i}.mask) = 0;
    glen = imgradient(input{i}.gx,input{i}.gy);
    glen_ex = imgradient(input{i}.gx_ex,input{i}.gy_ex);
    
    maxmag = max(glen(:));
    input{i}.gx = input{i}.gx ./ maxmag;
    input{i}.gy = input{i}.gy ./ maxmag;
    
    maxmag_ex = max(glen_ex(:));
    input{i}.gx_ex = input{i}.gx_ex ./ maxmag_ex;
    input{i}.gy_ex = input{i}.gy_ex ./ maxmag_ex;
        
    input{i}.grad(~input{i}.mask) = 0;
    maxmag = max(input{i}.grad(:));
    input{i}.grad = input{i}.grad ./ maxmag;
    
    maxmag_ex = max(input{i}.grad_ex(:));
    input{i}.grad_ex = input{i}.grad_ex ./ maxmag_ex;
    
    
    %---------fr_info------------
    p_z_c_mask = s_mask_contour(input{i}.mask);
    p_z_ex_only_mask = input{i}.mask_ex&~input{i}.mask;
    
    pval=double(input{i}.lab(s_ind2to3(p_z_c_mask)));
    pval_ex=double(input{i}.lab_ex(s_ind2to3(p_z_ex_only_mask)));
        
    pval_rgb=input{i}.rgb(s_ind2to3(p_z_c_mask));
    pval_ex_rgb=input{i}.rgb_ex(s_ind2to3(p_z_ex_only_mask));
    
    p_g_info=s_get_g_info(input{i},p_z_c_mask,p_z_ex_only_mask,options.nbins);

    fr_info{i}.c_mask = p_z_c_mask;
    fr_info{i}.ex_only_mask = p_z_ex_only_mask;
    fr_info{i}.c_vals = pval;
    fr_info{i}.ex_only_vals = pval_ex;
    fr_info{i}.c_vals_rgb = pval_rgb;
    fr_info{i}.ex_only_vals_rgb = pval_ex_rgb;
    fr_info{i}.g_info = p_g_info;
end

if options.is_use_palette_pxls_vals
    maxfrsz=max(cell2mat(cellfun(@(x) {x.origsz_ex},input)),[],1); 
    concat_rgb=(cell2mat(cellfun(@(x) {s_canvasSize(x.rgb_ex,maxfrsz)},input)));
    concat_lab=(cell2mat(cellfun(@(x) {s_canvasSize(x.lab_ex,maxfrsz)},input)));
    [labplt,map,rgbplt]=extractPalette(concat_rgb,options.palette_sz,concat_lab);
    for i=1:n_pics
        sti= (i-1)*maxfrsz(1)+1;
        eni=sti+maxfrsz(1)-1;
        curr_map = map(sti:eni,:);

        p_z_c_mask = s_mask_contour(input{i}.mask);
        [~,mask_i]=mask_indices(p_z_c_mask, maxfrsz);
        fr_info{i}.c_vals = labplt(curr_map(mask_i),:);
        fr_info{i}.c_vals_rgb = rgbplt(curr_map(mask_i),:);

        p_z_ex_only_mask = input{i}.mask_ex&~input{i}.mask;
        [~,mask_i]=mask_indices(p_z_ex_only_mask, maxfrsz);
        fr_info{i}.ex_only_vals = labplt(curr_map(mask_i),:);
        fr_info{i}.ex_only_vals_rgb = rgbplt(curr_map(mask_i),:);
    end
end
end

function newoptions = add_missing_options(oldoptions)
    newoptions=oldoptions;
    if ~isfield(oldoptions,'is_grad_dir_on_plt_output')
        newoptions.is_grad_dir_on_plt_output = false;
    end
    if ~isfield(oldoptions,'is_smooth_before_grad')
        newoptions.is_smooth_before_grad = true;
    end
    if ~isfield(oldoptions,'is_use_rand_frag_rot')
        newoptions.is_use_rand_frag_rot = false;
    end
    if ~isfield(oldoptions,'is_use_dist_colors_4grad')
        newoptions.is_use_dist_colors_4grad = true;
    end
    if ~isfield(oldoptions,'reduce_poly_length_to')
        newoptions.reduce_poly_length_to = 0.2;
    end
    if ~isfield(oldoptions,'nbins')
        newoptions.nbins = 80;
    end
    if ~isfield(oldoptions,'is_use_palette_pxls_vals')
        newoptions.is_use_palette_pxls_vals = false;
    end
    if ~isfield(oldoptions,'palette_sz')
        newoptions.palette_sz = 7;
    end
end