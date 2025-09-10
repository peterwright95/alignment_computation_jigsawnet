% function [result,execTimeMs]=assembly_main(pics_dir,seriesname,parts,ver,ext,rot_angles)
pics_dir='.'; % folder path that contains the fragments images
seriesname='hovav_x'; % name of the fragments series (<series_name>_p<frag_num><ver><ext>)
parts=[1:18]; % the fragment numbers range for the assembly calculation (e.g. 1:18)
ver='_m'; % the series ver (e.g. '_m')
ext='.png'; % the framgnet images extension (e.g. '.png')
rot_angles=round(rand(1,numel(parts))*359); % the fragments images initial random rotation, one for each fragment



try
result=[];
execTimeMs=[];

addpath(genpath(pwd));
tStart=tic;
%---------load fragments------------
if exist(pics_dir,'dir')
    addpath(pics_dir);
else
    error('pics directory not found');
end
% parts=str2double(parts);
full_pics = s_genNames(seriesname,parts,ver,ext);

%-----------init rotations----------
n_rot = 80;

dis_options.rot = linspace(0,360, n_rot + 1); dis_options.rot(end) = []; %remove 360 [93,207, 90,108];%
n_rot = numel(dis_options.rot); %reset in case rotations have been explicitly chosen
%-----------------------------------

load_options.is_grad_dir_on_plt_output = false;
load_options.is_smooth_before_grad = true;
load_options.is_use_rand_frag_rot = true; % rotate each fragment by a random angle
if numel(full_pics)~=numel(rot_angles)
    error('Not enough initial rotation angles supplied. expected=%d, recieved=%d',numel(full_pics),numel(rot_angles));
end
load_options.rand_angles = rot_angles;
load_options.is_use_dist_colors_4grad = true;
load_options.reduce_poly_length_to = 0.2;
load_options.nbins = n_rot; %number of angle bins for gradient histogram calculation
load_options.palette_sz = 10;
load_options.is_use_palette_pxls_vals = false;
[full_input, full_fr_info] = s_loadpics(full_pics, load_options);
%-----------------------------------
%-----init algorithm variables------
% out_folder = 'AssemblyRandRot'; %'ForPaper';'Assembly results';'AssemblyRandRot';'Classification'
% out_subfolder_mode = struct('mode','byseries','run_num',[]); %'bydate', 'byseries'
% [diary_base_filename,run_num]=outfilename( full_pics,n_rot,'txt',['assemb' '_log'],out_folder,out_subfolder_mode);
% out_subfolder_mode.run_num = run_num;
diary_base_filename = 'log/assembly_log.txt';
diary(ufilename(diary_base_filename));
fprintf('Loaded fragments: %s\n',string(join(full_pics,', ')));
fprintf('Fragments input rotations: %s\n',mat2str(cellfun(@(x) x.rand_angle,full_input)'));
dis_options.sampling_res = 5; % sampling resolution: number of pixels between each translation
dis_options.small_ov_percentage = 0.05; %controls threshold for dropping low count overlap transformations
dis_options.smoothing_deg_window = 13;

fprintf('----------------------\n');
fprintf('Init options: n_rot: %d, samp_res: %d\n',n_rot,dis_options.sampling_res);
fprintf('----------------------\n');
%% Dissimilarity
fprintf('-------------Starting dissimilarity--------------\n');
tStartDis=tic;
[full_dis_info]=s_calcerrors( full_input,full_fr_info,dis_options);
t=toc(tStartDis);fprintf('Total Dissimilarity Time: %d min, %0.3f sec\n',floor(t/60),mod(t,60));

%% Reduce number of tforms
dis_version=2;
fprintf('Using dissimilarity version: %d\n',dis_version);
n_reduced_tforms = 40000;
reduced_parts_i = 1:numel(full_pics);
input = full_input(reduced_parts_i);
fr_info = full_fr_info(reduced_parts_i);
[dis_info]=s_reduce_nsamples(n_reduced_tforms,reduced_parts_i,dis_version,full_dis_info,full_pics);

%% Pairwise Confidece
confoptions.iou_thresh = 0.85;%for top scores (intersection over union threshold)
confoptions.top_i=1;
confoptions.gamma_H=0.85;
confoptions.gamma_L=0.3;

fprintf('----------------------\n');
fprintf('Conf options: iou_thresh: %0.3f, gamma_L: %0.3f, gamma_H: %0.3f\n',confoptions.iou_thresh,confoptions.gamma_L,confoptions.gamma_H);
fprintf('----------------------\n');
fprintf('-------------Starting confidence--------------\n');
[ dis_info.conf,dis_info.best_match_info ] = confidence( dis_info,confoptions,input );

%% reassembly
sorted_1st_pair_i=1;
reassemb_options.iou_thresh = confoptions.iou_thresh;
reassemb_options.gamma_H = confoptions.gamma_H;
reassemb_options.gamma_L = confoptions.gamma_L;
reassemb_options.start_with_pair = sorted_1st_pair_i;
reassemb_options.max_output_images = 1;
reassemb_options.sampling_res = dis_options.sampling_res;
reassemb_options.dis_version = dis_version;
reassemb_options.small_ov_percentage = 0.1;
reassemb_options.smoothing_deg_window = dis_options.smoothing_deg_window;
reassemb_options.rot_refine_window=3;
reassemb_options.tr_refine_window=15;
reassemb_options.refine_alpha=0.0;
reassemb_options.reduce_poly_length_to=load_options.reduce_poly_length_to;
reassemb_options.save_top_conf_frag_num = [];
reassemb_options.is_output_assemb_progress = false;
reassemb_options.progress_out_folder=[];
fprintf('----------------------\n');
fprintf('Reassemb options: gamma_L: %0.3f, gamma_H: %0.3f, samp_res: %d, iou_thresh: %0.3f, sorted_1st_pair_i: %d\n'...
    ,reassemb_options.gamma_L,reassemb_options.gamma_H,reassemb_options.sampling_res,reassemb_options.iou_thresh,reassemb_options.start_with_pair);
fprintf('----------------------\n');
fprintf('-------------Starting reassembly--------------\n');
[fullrgb,finalTr,finalRot,order,totalconf,Adj,fullmask] = reassembly(input, fr_info, dis_info, reassemb_options);
t=toc(tStart);fprintf('Total Assembly Time: %d min, %0.3f sec\n',floor(t/60),mod(t,60));
execTimeMs = t;
fprintf('Final assembly transformations:\n');
disp([finalTr,finalRot]);
fprintf('Fragments assembly order (by fragment number): %s\n', mat2str(order{1}));
[neihbor_measure,correct_matches,total_matches,avgdistance,avg_iou]=d_neighbor_measure(full_pics,[finalTr,finalRot],cellfun(@(x) x.rand_angle,full_input)');
fprintf('Assemb properties: neighbor measure: %0.4f, (%d / %d), avg distance: %0.4f, avg iou: %0.4f, total conf: %0.4f\n', neihbor_measure, correct_matches, total_matches,avgdistance,avg_iou,totalconf(1));
diary('off');
result = fullrgb{1};
result_mask = fullmask{1};
result_after_rot_mask = imrotate(result_mask,-rot_angles(order{1}(1)));
result = imrotate(result,-rot_angles(order{1}(1)));
result(s_ind2to3(~result_after_rot_mask))=255;
result = s_imtrim(result, result_after_rot_mask);
% result = permute(result,[2,3,1]);
% for i=1:numel(fullrgb)
%     [~,outfn_wo_ext,outfn_ext]=fileparts(outfile);
%     uimwrite(fullrgb{i},fullfile(pics_dir,[outfn_wo_ext outfn_ext]));
% end
figure;imshow(result);
catch ME
    fprintf('Warning: An error occured, continue to next execution.\n\nException:\n%s\n',getReport(ME));
end
% end