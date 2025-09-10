function [ dis,dis_ov,n_ov,overlap_mask,ov_limits] = s_dissimilarity( fixed_colors,mov_colors,fixed_ex_only_i,mov_c_rc,fixed_c_i,tr,no_ov_val,sz,Mi,fixed_g_info,mov_g_info,small_ov_thresh,smoothing_deg_window)
alpha=2.0;
beta=2.0;
N=size(tr,1); % number of sample points
M=size(mov_c_rc,1); % number of contour points
T=numel(fixed_ex_only_i); % number of extrapolated only points

if ~isequal(size(fixed_g_info.ex_mag),[T,1]) ||...
   ~isequal(size(fixed_g_info.ex_bins),[T,1]) ||...
   ~isequal(size(mov_g_info.c_mag),[M,1]) ||...
   ~isequal(size(mov_g_info.c_bins),[M,1]) ||...
   mov_g_info.nbins ~= fixed_g_info.nbins ||...
   (~isempty(fixed_colors) && ~isequal(size(fixed_colors),[T,3])) ||...
   (~isempty(mov_colors) && ~isequal(size(mov_colors),[M,3]))
        error('input parameters size error');
end
NB=fixed_g_info.nbins;
tr = round(tr);

[shft_all_mov_cont_ind] = shift_points(mov_c_rc,tr,sz);%NxM, 1xM

[overlap_mask, overlap_loc] = ismember(shft_all_mov_cont_ind, fixed_ex_only_i); %NxM

n_ov = sum(overlap_mask,2); %Nx1

%remove physically impossible translations ( nonex overlapping)
if ~isempty(fixed_c_i)
    illegaloverlap_mask = any(ismember(shft_all_mov_cont_ind, fixed_c_i),2); %Nx1
else
    illegaloverlap_mask = false(N,1);%Nx1
end
tr2remove_mask = n_ov <= small_ov_thresh | illegaloverlap_mask; %Nx1
tr2calc_mask = ~tr2remove_mask;

dis_overlap = inf(N,1);
dis_v1 = inf(N,1);
dis_v2 = inf(N,1);
dis_v3 = inf(N,1);
ov_limits = nan(N,2);
    
if any(n_ov)
    tr2calc_ovmask = overlap_mask(tr2calc_mask,:);
    tr2calc_ovloc = overlap_loc(tr2calc_mask,:);
    [h_mov,h_fix]=calc_grad_hist(tr2calc_ovmask',tr2calc_ovloc',mov_g_info.c_bins,fixed_g_info.ex_bins,mov_g_info.c_mag,fixed_g_info.ex_mag,NB);
    h_mov=h_mov'+eps;
    h_fix=h_fix'+eps;
    h_mov=smoothdata(h_mov,2,'gaussian',round(smoothing_deg_window/mov_g_info.binwidth));
    mov_g_mag = sum(h_mov,2);
    h_fix=smoothdata(h_fix,2,'gaussian',round(smoothing_deg_window/fixed_g_info.binwidth));
    h_cityblock=abs(h_mov-h_fix); %N X NB
    h_jconv=h_cityblock.*log(max(h_mov,h_fix)./min(h_mov,h_fix)); %N X NB
    h_dist=sum(h_cityblock + h_jconv,2); %N X 1
    [overlap_mov_t,overlap_mov_i]=find(tr2calc_ovmask); % Kx1 : K = total overlapping points
    overlap_mov_t=tocolumn(overlap_mov_t);
    overlap_mov_i=tocolumn(overlap_mov_i);
    overlap_fixed_i=tocolumn(tr2calc_ovloc(tr2calc_ovmask));

    ov_dists = s_pdist(mov_colors(overlap_mov_i,:),fixed_colors(overlap_fixed_i,:),'1to1');
    tr2calc_dis_ov = accumarray(overlap_mov_t,ov_dists,[sum(tr2calc_mask),1],@sum,inf);
    tr2calc_n_ov = n_ov(tr2calc_mask);
    dis_overlap(tr2calc_mask) = tr2calc_dis_ov;
%     dis_v1(tr2calc_mask) = (tr2calc_dis_ov./tr2calc_n_ov) .* (1-tr2calc_n_ov./M).^alpha;
%     dis_v2(tr2calc_mask) = ((tr2calc_dis_ov+h_dist)./tr2calc_n_ov) .* (1-tr2calc_n_ov./M).^alpha;
    dis_v1(tr2calc_mask) = (tr2calc_dis_ov) ./ (tr2calc_n_ov.^alpha);
    dis_v2(tr2calc_mask) = (tr2calc_dis_ov+h_dist) ./ (tr2calc_n_ov.^alpha);
    dis_v3(tr2calc_mask) = (tr2calc_dis_ov+h_dist) ./ (tr2calc_n_ov.^alpha) + 1./(mov_g_mag.^beta);
%     dis_v3(tr2calc_mask) = (tr2calc_dis_ov+h_dist) ./ (beta.*tr2calc_n_ov+(1-beta).*mov_g_mag).^alpha;
    ov_limits(tr2calc_mask,:)=calc_ov_limits(tr2calc_ovmask',Mi)';
%     dis_v3 = dis_overlap;
%     dis_v3 = dis_v1+dis_v2;%(dis_overlap+h_dist2) ./ (n_ov.^alpha);
%     dis_v3 = (dis_overlap + no_ov_val*(M-n_ov));
    
end


dis_v1(tr2remove_mask) = inf;
dis_v2(tr2remove_mask) = inf;
dis_v3(tr2remove_mask) = inf;
dis_overlap(tr2remove_mask) = inf;
n_ov(tr2remove_mask) = 0;
dis_ov(:,1)=dis_overlap;
dis_ov(:,2)=dis_overlap;
dis_ov(:,3)=dis_overlap;
dis(:,1)=dis_v1;
dis(:,2)=dis_v2;
dis(:,3)=dis_v3;
end
function [shft_ind,shft_alltr_points_r,shft_alltr_points_c] = shift_points(points_rc,tr,sz)
N=size(tr,1); % number of sample points
M=size(points_rc,1); % number of contour points
% alltr_points_rc = repmat(permute(points_rc,[3,1,2]),N,1,1); %NxMx2
% shft_alltr_points_rc = round(alltr_points_rc + repmat(permute(tr,[1,3,2]),1,M,1)); %NxMx2
% shft_alltr_points_rc = round(cat(3,bsxfun(@plus,tr(:,1),points_rc(:,1)'),bsxfun(@plus,tr(:,2),points_rc(:,2)'))); %NxMx2
shft_alltr_points_rc = round(bsxfun(@plus,permute(tr,[1,3,2]),permute(points_rc,[3,1,2])));%NxMx2
% if any(any(shft_alltr_points_rc(:,:,1)<1)) || any(any(shft_alltr_points_rc(:,:,1)>sz(1))) ||...
%     any(any(shft_alltr_points_rc(:,:,2)<1)) || any(any(shft_alltr_points_rc(:,:,2)>sz(2)))
%     error('shifted indices out of range');
% end
shft_alltr_points_r = shft_alltr_points_rc(:,:,1);%NxM
shft_alltr_points_c = shft_alltr_points_rc(:,:,2);%NxM
illigalpoints_mask = shft_alltr_points_r<1 | shft_alltr_points_r>sz(1) | shft_alltr_points_c<1 | shft_alltr_points_c>sz(2);%NxM
if any(illigalpoints_mask(:))
    warning('%d/%d of the shifted indices are out of range : [1:%d,1:%d].',nnz(illigalpoints_mask(:)),N*M,sz(1),sz(2));
end
shft_alltr_points_r(illigalpoints_mask)=1;
shft_alltr_points_c(illigalpoints_mask)=1;
shft_ind = sub2ind(sz(1:2),shft_alltr_points_r,shft_alltr_points_c); %NxM
shft_ind(illigalpoints_mask)=0;
% non_shft_ind = sub2ind(sz(1:2),points_rc(:,1),points_rc(:,2))'; %1xM
end
function [weights] = grad_weights(grad_w8s)
if iscolumn(grad_w8s)
    grad_w8s=grad_w8s';
end
epsv=0.01;
grad_vals = grad_w8s + epsv; %1xM
grad_vals = grad_vals ./ sum(grad_vals);
weights = grad_vals;%1xM
end
function debug_showoverlap(tmpi,nfig,shft_all_mov_cont_ind,fix_ex_ind,overlap_mask,tr,sz)
figure(nfig);
tmp=zeros(sz(1:2));
tmp(shft_all_mov_cont_ind(tmpi,:)) = 0.5;
tmp(fix_ex_ind) =tmp(fix_ex_ind) +0.3;
imshow(tmp);
fprintf('overlap: %d, shift: [%d,%d]\n',sum(overlap_mask(tmpi,:),2),tr(tmpi,1),tr(tmpi,2));
end