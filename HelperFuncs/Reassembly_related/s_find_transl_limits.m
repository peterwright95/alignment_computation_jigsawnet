function [ limits,cn_offset ] = s_find_transl_limits( f, T )
NV = size(T,1);
if ~iscell(f) || numel(f)~= NV
    error('There should be translation for each fragment, numel(f)=%d,number of T=%d',numel(f),NV);
end
% minloc = [inf,inf];
% maxloc = [-inf,-inf];
% for v = 1:NV
% cent = size(f{v}.mask)/2;
% a = [T(v,1:2);T(v,1:2)] + [1,1;size(f{v}.mask)] - [cent;cent];
% minloc = min(minloc,min(a,[],1));
% maxloc = max(maxloc,max(a,[],1));
% end
% limits = ceil([2*max(abs([maxloc(1),minloc(1)])),2*max(abs([maxloc(2),minloc(2)]))]);
Tbbox_loc = zeros(NV,4);
for v = 1:NV
    after_rot_bbox_loc = ceil(abs(rot2D(T(v,3)))*size(f{v}.mask)')';%1x2
    voffset = after_rot_bbox_loc./2;
    Tbbox_loc(v,:) = [1,1,after_rot_bbox_loc]-[voffset,voffset] + [T(v,1:2),T(v,1:2)];
end
llimits = floor(min(Tbbox_loc(:,1:2),[],1));
ulimits = ceil(max(Tbbox_loc(:,3:4),[],1));
limits_offset = abs(llimits) + 1;
limits = ulimits + limits_offset;
cn_offset = -(ulimits + llimits)/2;
end

