function outInfo = generate_ext_for_folder(folderPath)
% simple_extrapolate_folder  Minimal wrapper to create _m and _m_ext images.
%
% Key difference from earlier draft: we NO LONGER rescale fragments.
%   resize_factor is forced to 1 so the produced <name>_m.png and
%   <name>_m_ext.png share the same pixel scale as the original fragment.
%   (They will still have the same padded canvas internally added by the
%   original extrapolation code; _m and _m_ext have identical size, _m_ext
%   just fills an extrapolated halo.)
%
% Why: User already has final cropped fragments and does not want any
% geometric distortion introduced by normalization. Setting resize=1 keeps
% fragment shapes unchanged (aside from the uniform padding the original
% pipeline always adds around both outputs).
%
% This mimics the repository behavior where a fragment and its _ext share
% the same dimensions. If a <name>_full.png twin does not exist (the
% original pipeline expected one), we create a duplicate so that
% extern_imextrapolate can read <name>_full.png without errors.
%
% Usage:
%   out = simple_extrapolate_folder('Pictures/MIT_ex');
%
% Inputs:
%   folderPath : directory containing fragment images (*.png). Each base
%                image will produce <name>_m.png and <name>_m_ext.png in
%                the SAME directory (no subfolder).
%
% Output (struct outInfo):
%   .resize_factor   - always 1
%   .mode            - 'm'
%   .input_images    - cell array of processed base image filenames
%   .skipped         - images skipped because _m_ext already existed
%   .created_full    - list of *_full duplicates we generated
%   .out_dirs        - output directories returned by extrapolate_fragments
%   .out_img_postfix - postfix returned (e.g., '_m')
%
% Notes:
%   - Single mandatory argument (folderPath).
%   - Delete existing *_m_ext.* to force regeneration.
%   - Only *.png supported for simplicity.
%   - If you want to eliminate the internal padding too, you'd need to
%     modify extern_imextrapolate (it pads via s_canvasSize). Ask if wanted.
%
if ~isfolder(folderPath)
    error('Folder not found: %s',folderPath);
end

pattern = '*.png';
mode = 'm';  % kept for naming (_m / _m_ext) but NOT used as subfolder now
ext_postfix = '_ext';

files = dir(fullfile(folderPath, pattern));
if isempty(files)
    warning('No PNG files found in %s', folderPath);
end

% Filter base images (exclude ones already ending with _m.png or _m_ext.png)
base = {};
skipped = {};
for i=1:numel(files)
    [~,nm,ext] = fileparts(files(i).name);
    if endsWith(nm, ['_' mode]) || endsWith(nm, ['_' mode ext_postfix])
        continue; % already derived artifact
    end
    % Check if extended exists already
    extExists = exist(fullfile(folderPath, [nm '_' mode '_ext' ext]), 'file');
    if extExists
        skipped{end+1} = files(i).name; %#ok<AGROW>
        continue;
    end
    base{end+1} = files(i).name; %#ok<AGROW>
end

if isempty(base)
    fprintf('Nothing to process. All images appear to be already extrapolated or none found.\n');
    outInfo = struct('resize_factor',[], 'mode',mode, 'input_images',{base}, 'skipped',{skipped}, 'out_dirs',{{}}, 'out_img_postfix','');
    return;
end

% Force no geometric rescaling.
resize_factor = 1;

% Ensure expected *_full.png exists (duplicate original if missing).
created_full = {};
for i=1:numel(base)
    [~,nm,ext] = fileparts(base{i}); %#ok<ASGLU>
    fullName = fullfile(folderPath, [nm '_full' ext]);
    if ~exist(fullName,'file')
        copyfile(fullfile(folderPath, base{i}), fullName);
        created_full{end+1} = [nm '_full' ext]; %#ok<AGROW>
    end
end

% Add folder so extrapolate_fragments can find images by base name
addpath(folderPath);

% Provide pic names (without path) so *_full logic works via MATLAB path.
pics = base; %#ok<NASGU>

% Call existing extrapolation pipeline
[out_dirs, out_img_postfix] = extrapolate_fragments(pics, resize_factor, mode);

% Move generated outputs (placed by extrapolate_fragments in a 'm' subdir)
% up into the original folderPath.
produced_dir_primary = out_dirs{1};
for i=1:numel(base)
    [~,nm,ext] = fileparts(base{i});
    src_main = fullfile(produced_dir_primary, [nm out_img_postfix ext]);
    src_ext  = fullfile(produced_dir_primary, [nm out_img_postfix '_ext' ext]);
    dst_main = fullfile(folderPath, [nm out_img_postfix ext]);
    dst_ext  = fullfile(folderPath, [nm out_img_postfix '_ext' ext]);
    if exist(src_main,'file'); movefile(src_main, dst_main, 'f'); end
    if exist(src_ext,'file');  movefile(src_ext,  dst_ext,  'f'); end
end
% Attempt to remove now-empty produced_dir_primary and its parent if empty.
try
    dlist = dir(produced_dir_primary);
    if numel(dlist) <= 2
        rmdir(produced_dir_primary);
    end
catch, end %#ok<CTCH>
out_dirs = {folderPath};

outInfo = struct('resize_factor', resize_factor, ...
              'mode', mode, ...
              'input_images', {base}, ...
              'skipped', {skipped}, ...
              'created_full', {created_full}, ...
              'out_dirs', {out_dirs}, ...
              'out_img_postfix', out_img_postfix);

fprintf(['Extrapolation complete. resize_factor=%0.2f | processed=%d | ' ...
        'skipped(existing)=%d | created_full=%d\n'], ...
    resize_factor, numel(base), numel(skipped), numel(created_full));
end


