clear variables;
masks_folder = 'pictures/Masks/s8';
masks_ext = '.jpg';
orig_img_filenames = {'input.jpg'};
out_series_names = {'hovav_x'};
out_fileext = '.png';

addpath(genpath(pwd));
for i=1:numel(out_series_names)
out_series_name=out_series_names{i};
orig_img_filename=orig_img_filenames{i};
[fr_out_dir,n_fragments,avg_img_size]=img2fragments(masks_folder,masks_ext,orig_img_filename,out_series_name,out_fileext);
path(fr_out_dir,path);

%%
desired_img_size=[250 250];
resize_factor = min(1,round(mean(desired_img_size./avg_img_size.*100))./100);
extrap_mode = 'm';
parts=1:n_fragments;
extrap_pics=s_genNames(out_series_name,parts,'',out_fileext);

[extrap_fr_outdirs,out_img_postfix]=extrapolate_fragments(extrap_pics,resize_factor,extrap_mode);
path(extrap_fr_outdirs{1},path);

outpics=s_genNames(out_series_name,parts,out_img_postfix,out_fileext);
gtT=d_gt_for_masks(masks_folder,orig_img_filename,resize_factor);
str=d_Ttojson(gtT,outpics);
fileID = fopen(fullfile(extrap_fr_outdirs{1},['reduced_size_by_' num2str(resize_factor) '.txt']),'w');
fprintf(fileID,str);
fclose(fileID);
gt_img=d_Tforms2im(gtT,s_loadpics(outpics));
uimwrite(gt_img,fullfile(fileparts(fr_out_dir),[out_series_name,out_img_postfix,'_gt.jpg']));
end
% gt_json_filename = 'groundtruth.json';
% fileID = fopen(gt_json_filename,'rt');
% gt_st=jsondecode(fscanf(fileID,'%s'));
% fclose(fileID);
% fileID = fopen(which(gt_json_filename),'wt');
% new_example_gt_st=jsondecode(str);
% gt_st(end+1)=new_example_gt_st;
% fprintf(fileID,jsonencode(gt_st));
% fclose(fileID);