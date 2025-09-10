function [ G, Gx, Gy ] = s_palette_grad( rgb,lab, palette_sizes,method,is_use_dist_colors,is_smooth_before_grad )
%S_PALETTE_GRAD gradient on image in lab space with palette extraction
sigma = 1;
if ~exist('palette_sizes','var')
    palette_sizes = [3,4,6,8,10];
end
if ~exist('method','var')
    method = 'sobel';
end
if ~exist('is_use_dist_colors','var')
    is_use_dist_colors = false;
end
if ~exist('is_smooth_before_grad','var')
    is_smooth_before_grad = false;
end
if (is_smooth_before_grad)
   rgb = imgaussfilt(rgb,sigma);
   lab = imgaussfilt(lab,sigma);
end
G = zeros(size(rgb,1),size(rgb,2));
if nargout > 1
    Gx = zeros(size(rgb,1),size(rgb,2));
    Gy = zeros(size(rgb,1),size(rgb,2));
end
for psz=palette_sizes
    [labplt,map,rgbplt]=extractPalette(rgb,psz,lab);
    if is_use_dist_colors
        bgc = get_background_color(reshape(labplt(map(:), :),size(map,1), size(map,2),3));
        dists=pdist2(bgc,labplt);
        bgi = find(dists==0);
        colors = zeros(size(labplt));
        colors(bgi,:) = labplt(bgi,:);
        [~,colors(setdiff(1:psz,bgi),:)] = dist_colors(psz-1,rgbplt(bgi,:)/255,@rgb2lab);
    else
        colors = labplt;
    end
    pal_lab = reshape(colors(map(:), :),size(map,1), size(map,2),3);
    assert(isa(pal_lab,'double'),'palette expected to be double');
    if nargout>1
        [Ix,Iy]=imgradientxy(pal_lab(:,:,1),method);
        Gx = Gx + Ix./psz;
        Gy = Gy + Iy./psz;
        I=imgradient(Ix,Iy);
    else
        I=imgradient(pal_lab(:,:,1),method);
    end
    G=G+I./psz;
end


end

