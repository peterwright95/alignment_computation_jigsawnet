function [outdirs,out_img_postfix] = extrapolate_fragments(pics,resize_factor,mode)

% clear variables;close all;
% pics=[s_genNames('nicholas_s6',1:42,'','.png')];
% resize_factor = 1.0;
load_img_postfix='_full';
% mode='m';
out_img_postfix = ['_' mode];
im2bin_mode='other'; %other/palette
outdirs = cell(numel(pics),1);
for j=1:1
    ex = cell(numel(pics),1);
    img = cell(numel(pics),1);
    targetregion = cell(numel(pics),1);
    imgCompleteFinal = cell(numel(pics),1);
    erosion=j;
    for i=1:numel(pics)
        [ ex{i},img{i},targetregion{i},imgCompleteFinal{i} ] = extern_imextrapolate(pics{i},mode,load_img_postfix,resize_factor,im2bin_mode,erosion);

        [~,picname,picext]=fileparts(pics{i});
        [dir,~,~]=fileparts(which([picname,load_img_postfix,picext]));
        if ~strcmp(load_img_postfix,'')
            outdir = fileparts(dir);
        else
            outdir = dir;
        end
%         uimwrite(ex{i}, fullfile(outdir,mode,['erosion_' num2str(erosion)], [picname,out_img_postfix,['e' num2str(erosion)], '_ext', picext]));
%         uimwrite(img{i}, fullfile(outdir,mode,['erosion_' num2str(erosion)], [picname,out_img_postfix,['e' num2str(erosion)], picext]));
        uimwrite(ex{i}, fullfile(outdir,mode,[picname,out_img_postfix, '_ext', picext]));
        uimwrite(img{i}, fullfile(outdir,mode,[picname,out_img_postfix, picext]));
        outdirs{i}=fullfile(outdir,mode);
        % figure(i);
        % subplot(1,4,1);imshow(img{i});title('original');
        % subplot(1,4,2);imshow(imgCompleteFinal{i});title('whole image completion');
        % ex_with_cont = ex{i};
        % ex_with_cont(s_ind2to3(s_mask_contour(targetregion{i}))) = 0;
        % subplot(1,4,3);imshow(ex_with_cont);title('extrapolated with contour');
        % subplot(1,4,4);imshow(ex{i});title('extrapolated');
    end
end
outdirs=unique(outdirs);

