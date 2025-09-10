function [ fullrgb,fullmask ] = d_Tforms2im( T,f,isTrim,is_show_c_color,c_size,c_uniform_color)
%S_TFORMS2IM construct the final image from transformations
%   T: NVx3, list of all transformations
%   f: NVx1 cell array, with all the input information
if ~exist('isTrim','var')
    isTrim = true;
end
if ~exist('c_uniform_color','var')
    c_uniform_color = [];
end
if ~exist('c_size','var')
    c_size = 2;
end
if ~exist('is_show_c_color','var')
    is_show_c_color = false;
end

NV = size(T,1);
[finalsz,cn_offset] = s_find_transl_limits( f, T ); %expand limits for rotation
fullrgb = zeros([finalsz,3],'uint8');
fullmask = false(finalsz);
T(:,1:2) = T(:,1:2)+cn_offset;

shifted_masks = cell(NV,1);
for v = 1:NV
    [fullrgb,fullmask,~,shifted_masks{v}] = s_imcomb(fullrgb,f{v}.rgb,fullmask,f{v}.mask,T(v,:),'assrc');
end

if is_show_c_color
    if ~isempty(c_uniform_color)
        vcolors=repmat(c_uniform_color,NV,1);
    else
        [~,~,assemb_colors_plt]=extractPalette(fullrgb,7);
        assemb_colors_plt = double(assemb_colors_plt)/255;
        vcolors = dist_colors(NV,assemb_colors_plt)*255;
    end

    for v = 1:NV
        vcolor = vcolors(v,:);
        vc_mask = shifted_masks{v} & ~s_scale_mask(shifted_masks{v}, -c_size);
        vc_mask_i = s_ind2to3(vc_mask);
        fullrgb(vc_mask_i) = repmat(vcolor,size(vc_mask_i,1),1);
    end
end
if (isTrim)
    [fullrgb,~,fullmask] = s_imtrim(fullrgb,fullmask);
end
end

