function [outdir,n_masks,avg_img_size]=img2fragments(masks_folder,masks_ext,orig_img_filename,out_series_name,out_fileext)

% pics = dir('C:\Users\nderech\Google Drive\Technion\Pictures\DryMud_fragments\s6\*.jpg');
% pics = dir('C:\Users\NivD\Google Drive\Technion\Pictures\DryMud_fragments\s6\*.jpg');
% pics = dir('D:\Projects\TempProjects\assemb\Pictures\DryMud_fragments\s9\*.jpg');
pics = dir(fullfile(masks_folder,['*',masks_ext]));

% pics = dir('C:\Users\nderech\Google Drive\Technion\Pictures\Frescoes\romanvilla\masks\*.jpg');
masks = arrayfun(@(x) {imbinarize(rgb2gray(imread(x.name)))},pics);
n_masks = numel(masks);
masks_sizes = cell2mat(cellfun(@(x) {size(x)},masks'));
if any(masks_sizes(:,1)~= masks_sizes(1,1)) || any(masks_sizes(:,2)~= masks_sizes(1,2))
    error('masks have different sizes. sizes: %s',mat2str(masks_sizes));
end
masks_sz = masks_sizes(1,:);

img_filename = orig_img_filename;% 'st_nicholas_roof.jpg';
img = imread(img_filename);
img_lab = rgb2lab(img);
img_sz = size(img);

if ~isequal(img_sz(1:2),masks_sz)
    if (img_sz(1)>img_sz(2) && masks_sz(1)<masks_sz(2)) ||...
       (img_sz(1)<img_sz(2) && masks_sz(1)>masks_sz(2))
        masks = cellfun(@(x) {imresize(x',img_sz(1:2))},masks);    
    else
        masks = cellfun(@(x) {imresize(x,img_sz(1:2))},masks);
    end    
end

unionmask = false(img_sz(1:2));
for i=1:n_masks
    unionmask=unionmask|masks{i};
end

series = out_series_name;
outdir = fullfile(fileparts(which(img_filename)),series,'full');
avg_img_size=zeros(1,2);
for i=1:n_masks
    wfrag = img;
    wfrag_lab = img_lab;
    bfrag = img;
    bfrag_lab = img_lab;
    outside = s_ind2to3(~masks{i});
    bfrag(outside) = zeros(size(outside));
    wfrag(outside) = zeros(size(outside))+255;
    bfrag_lab(outside) = zeros(size(outside));
    wfrag_lab(outside) = repmat([100,0,0],size(outside,1),1);
    mask_from_black_bg = s_im2bin(bfrag,bfrag_lab,'other');
    mask_from_white_bg = s_im2bin(wfrag,wfrag_lab,'other');
    bbg_iou=nnz(mask_from_black_bg&masks{i})/nnz(mask_from_black_bg|masks{i});
    wbg_iou=nnz(mask_from_white_bg&masks{i})/nnz(mask_from_white_bg|masks{i});
    if bbg_iou>wbg_iou
        tfrag = s_imtrim(bfrag,masks{i});
    else
        tfrag = s_imtrim(wfrag,masks{i});
    end
    avg_img_size=avg_img_size+size(tfrag(:,:,1));
    outfile = fullfile(outdir,[series,'_p',num2str(i),'_full',out_fileext]);
    uimwrite(tfrag,outfile);
end
avg_img_size=avg_img_size/n_masks;