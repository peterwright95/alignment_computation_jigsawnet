function [ gtT] = d_gt_for_masks( masks_folder, orig_img_file, orig_resize_factor)
%D_GT_FOR_MASKS Summary of this function goes here
%   Detailed explanation goes here

pics = dir(fullfile(masks_folder,'*.jpg'));
masks = arrayfun(@(x) {imbinarize(rgb2gray(imread(x.name)))},pics);

masks_sizes = cell2mat(cellfun(@(x) {size(x)},masks'));
if any(masks_sizes(:,1)~= masks_sizes(1,1)) || any(masks_sizes(:,2)~= masks_sizes(1,2))
    error('masks have different sizes. sizes: %s',mat2str(masks_sizes));
end
masks_sz = masks_sizes(1,:);

img_filename = orig_img_file;
img = imread(img_filename);
img_sz = size(img);

if ~isequal(img_sz(1:2),masks_sz)
    if (img_sz(1)>img_sz(2) && masks_sz(1)<masks_sz(2)) ||...
       (img_sz(1)<img_sz(2) && masks_sz(1)>masks_sz(2))
        masks = cellfun(@(x) {imresize(x',round(img_sz(1:2)*orig_resize_factor))},masks);    
    else
        masks = cellfun(@(x) {imresize(x,round(img_sz(1:2)*orig_resize_factor))},masks);
    end    
end
unionmask = false(size(masks{1}));
for i=1:numel(masks)
    unionmask=unionmask|masks{i};
end

gtT=zeros(numel(pics),3);
bb1=get_bounding_box(masks{1});
cn1=[bb1.rmin+bb1.rmax,bb1.cmin+bb1.cmax]/2;
for i=2:numel(pics)
    bb=get_bounding_box(masks{i});
    cn=[bb.rmin+bb.rmax,bb.cmin+bb.cmax]/2;
    gtT(i,1:2)=cn-cn1;
end
end

