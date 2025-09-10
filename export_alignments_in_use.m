pics_dir   = fullfile('Pictures','MIT_ex');
seriesname = 'fragment';
parts      = 1:9;
ver        = '_m';
ext='.png'; % image extension
% NOTE: Input images already have desired orientations; disable random rotation.
rot_angles=[]; % unused (kept for compatibility if re-enabled)
output_txt='alignments.txt'; % output file for pairwise alignments
top_k = 3; % number of transforms per pair to export
save_pair_images = true; % save combined images for each exported transform
pair_img_dir = fullfile('Pictures','PairAlignments');

try
addpath(genpath(pwd));
if exist(pics_dir,'dir')
    addpath(pics_dir);
else
    error('pics directory not found');
end
full_pics = s_genNames(seriesname,parts,ver,ext);

% --- init rotations (exactly as assembly_main) ---
n_rot = 80;
dis_options.rot = linspace(0,360,n_rot+1); dis_options.rot(end)=[]; % remove 360
n_rot = numel(dis_options.rot);

% --- load fragments (same fields/order) ---
load_options.is_grad_dir_on_plt_output = false;
load_options.is_smooth_before_grad = true;
load_options.is_use_rand_frag_rot = false; % DO NOT apply extra random rotations
load_options.is_use_dist_colors_4grad = true;
load_options.reduce_poly_length_to = 0.2;
load_options.nbins = n_rot;
load_options.palette_sz = 10;
load_options.is_use_palette_pxls_vals = false;
[full_input, full_fr_info] = s_loadpics(full_pics, load_options);

diary_base_filename = 'log/assembly_log.txt';
diary(ufilename(diary_base_filename));
fprintf('Loaded fragments: %s\n',string(join(full_pics,', ')));
% Rotations all zeros because we disabled random rotation
fprintf('Fragments input rotations disabled (all zero).\n');
dis_options.sampling_res = 5;
dis_options.small_ov_percentage = 0.05;
dis_options.smoothing_deg_window = 13;
fprintf('----------------------\n');
fprintf('Init options: n_rot: %d, samp_res: %d\n',n_rot,dis_options.sampling_res);
fprintf('----------------------\n');

% --- dissimilarity ---
fprintf('-------------Starting dissimilarity--------------\n');
full_dis_info = s_calcerrors(full_input,full_fr_info,dis_options);

% --- reduce number of tforms (same defaults) ---
dis_version=2; fprintf('Using dissimilarity version: %d\n',dis_version);
n_reduced_tforms = 40000;
reduced_parts_i = 1:numel(full_pics);
input = full_input(reduced_parts_i); %#ok<NASGU>
fr_info = full_fr_info(reduced_parts_i); %#ok<NASGU>
dis_info = s_reduce_nsamples(n_reduced_tforms,reduced_parts_i,dis_version,full_dis_info,full_pics);

% --- confidence ---
confoptions.iou_thresh = 0.85;
confoptions.top_i=1;
confoptions.gamma_H=0.85;
confoptions.gamma_L=0.3;
fprintf('----------------------\n');
fprintf('Conf options: iou_thresh: %0.3f, gamma_L: %0.3f, gamma_H: %0.3f\n',confoptions.iou_thresh,confoptions.gamma_L,confoptions.gamma_H);
fprintf('----------------------\n');
fprintf('-------------Starting confidence--------------\n');
[dis_info.conf,dis_info.best_match_info,~,allmatch_info] = confidence(dis_info,confoptions,full_input);

% --- export ONE best transform per pair (mirrors internal data, no extra ranking logic) ---
fid=fopen(output_txt,'w'); if fid<0, error('Cannot open %s',output_txt); end
cleanupObj=onCleanup(@() fclose(fid)); %#ok<NASGU>
% Node headers (zero-based)
for i=1:numel(full_pics); fprintf(fid,'Node %d\n',i-1); end
fprintf('Writing top-%d pairwise transforms to %s\n',top_k,output_txt);
if save_pair_images && ~exist(pair_img_dir,'dir'); mkdir(pair_img_dir); end
for pi = 1:size(dis_info.pairs,1)
    if pi>numel(allmatch_info) || isempty(allmatch_info{pi}) || ~isfield(allmatch_info{pi},'tform_i1')
        continue;
    end
    tlist = allmatch_info{pi}.tform_i1;
    if isempty(tlist); continue; end
    p = dis_info.pairs(pi,:); p1 = p(1); p2 = p(2);
    k = min(top_k, numel(tlist));

    for ki = 1:k
        ti = tlist(ki);
        if ti<1 || ti>size(dis_info.tforms{pi},1), continue; end

        % --- pose from confidence (dst p2 -> src p1) ---
        T   = dis_info.tforms{pi}(ti,:);   % [dy dx rdeg]
        dy  = T(1);  dx = T(2);  th = deg2rad(T(3));
        c   = cos(th); s = sin(th);

        % row/col rigid pieces (Y=row, X=col)
        Rrc = [ c  s  0;   -s  c  0;   0 0 1 ];       % rotation
        Trc = [ 1  0  dy;   0  1  dx;  0 0 1 ];       % translation

        % rotate about the center of dst (p2)
        [h2,w2,~] = size(full_input{p2}.rgb);
        cy = (h2-1)/2;  cx = (w2-1)/2;                % center in row/col coords
        Cto   = [1 0 cy; 0 1 cx; 0 0 1];
        Cfrom = [1 0 -cy; 0 1 -cx; 0 0 1];

        A = Trc * Cto * Rrc * Cfrom;                  % final 3×3, dst->src (row/col)

        ov_cnt = min(dis_info.n_ov1{pi}(ti), dis_info.n_ov2{pi}(ti));

        % --- seam coords (keep as before; src-frame points) ---
        try
            [coords_x,coords_y] = local_overlap_coords(full_input{p1}.mask, full_input{p2}.mask, dx, dy, T(3));
        catch
            coords_x = []; coords_y = [];
        end

        % --- emit line: IDs, ov_cnt, 9 elems of A (row-major), then 'line' + coords ---
        fprintf(fid, '%d %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f line', ...
            p1-1, p2-1, ov_cnt, ...
            A(1,1), A(1,2), A(1,3), ...
            A(2,1), A(2,2), A(2,3), ...
            A(3,1), A(3,2), A(3,3));

        for ci = 1:numel(coords_x)
            fprintf(fid, ' %.6f %.6f', coords_x(ci), coords_y(ci));  % x y == col row (src frame)
        end
        fprintf(fid, '\n');

        % --- optional preview image (unchanged; s_imcomb expects [dy dx rdeg]) ---
        if save_pair_images
            try
                comb = s_imcomb(full_input{p1}.rgb, full_input{p2}.rgb, ...
                                full_input{p1}.mask, full_input{p2}.mask, T);
                outfn = sprintf('pair_%02d_%02d_k%d.png', p1, p2, ki);
                imwrite(comb, fullfile(pair_img_dir, outfn));
            catch imgME
                fprintf('Warn: failed saving image for pair (%d,%d) k=%d: %s\n', p1, p2, ki, imgME.message);
            end
        end
    end
end
fprintf('Export complete.\n');
diary('off');
catch ME
    fprintf('Error: %s\n',ME.message);
end

function [coords_x,coords_y] = local_overlap_coords(mask1,mask2,dx,dy,rdeg)
% Rotate second mask by rdeg, shift by (dx,dy) relative to mask1, return overlap pixel coords
rot_mask2 = imrotate(mask2,rdeg);
finalsz = max(size(mask1),size(rot_mask2));
mask1c = s_canvasSize(mask1,finalsz);
rot_mask2c = s_canvasSize(rot_mask2,finalsz);
% shift coordinates of mask2 pixels
[r2,c2] = find(rot_mask2c);
shifted_r = r2 + dy; shifted_c = c2 + dx;
valid = shifted_r>=1 & shifted_r<=finalsz(1) & shifted_c>=1 & shifted_c<=finalsz(2);
shifted_r = shifted_r(valid); shifted_c = shifted_c(valid);
lin = sub2ind(finalsz,shifted_r,shifted_c);
shifted_mask2 = false(finalsz); shifted_mask2(lin) = true;
overlap = mask1c & shifted_mask2;
[coords_y,coords_x] = find(overlap);
% subsample if too many points
MAX_COORDS = 400;
if numel(coords_x) > MAX_COORDS
    step = ceil(numel(coords_x)/MAX_COORDS);
    coords_x = coords_x(1:step:end);
    coords_y = coords_y(1:step:end);
end
coords_x = double(coords_x); coords_y = double(coords_y);
end
