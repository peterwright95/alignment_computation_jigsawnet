% function [ req_pairs_i,dis1,dis2,gsim1,gsim2,tforms,n_ov1,n_ov2] =...
%     s_reduce_nsamples( req_n_reduced_tforms,dis_version,pairs,reduced_parts_idx,full_dis1,full_dis2,full_gsim1,full_gsim2,full_tforms,full_n_ov1,full_n_ov2 )
function [dis_info_out]=s_reduce_nsamples(req_n_reduced_tforms,reduced_parts_idx,dis_version,dis_info_in,full_pics)
[req_pairs_i,~] = find(all(ismember(dis_info_in.pairs,reduced_parts_idx),2));
min_n_tforms=min(cellfun(@(x) numel(x(:,dis_version)) ,dis_info_in.dis1(req_pairs_i)));
n_reduced_tforms = min(req_n_reduced_tforms,min_n_tforms);

srtdis_i = cellfun(@(dis1,dis2) {sorted_ind(dis1(:,dis_version)+dis2(:,dis_version),1,'ascend')},dis_info_in.dis1(req_pairs_i),dis_info_in.dis2(req_pairs_i));
% srtdis_i = cellfun(@(dis1,dis2) {sorted_ind(dis1(:,dis_version),1,'ascend')},dis_info_in.dis1(req_pairs_i),dis_info_in.dis2(req_pairs_i));
% dis_rti=cellfun(@(srt_dis_pi) {sort(srt_dis_pi(1:n_reduced_tforms,:),1,'ascend')},srtdis_i);
dis_rti=cellfun(@(srt_dis_pi) {(1:numel(srt_dis_pi))'},srtdis_i);
% dis_info_out.dis_ov1 = cellfun(@(x,f) {f(x,dis_version)},dis_rti,dis_info_in.dis_ov1(req_pairs_i));
% dis_info_out.dis_ov2 = cellfun(@(x,f) {f(x,dis_version)},dis_rti,dis_info_in.dis_ov2(req_pairs_i));
dis_info_out.dis1 = cellfun(@(x,f) {f(x,dis_version)},dis_rti,dis_info_in.dis1(req_pairs_i));
dis_info_out.dis2 = cellfun(@(x,f) {f(x,dis_version)},dis_rti,dis_info_in.dis2(req_pairs_i));
dis_info_out.tforms = cellfun(@(x,f) {f(x,:)},dis_rti,dis_info_in.tforms(req_pairs_i));
dis_info_out.n_ov1 = cellfun(@(x,f) {f(x)},dis_rti,dis_info_in.n_ov1(req_pairs_i));
dis_info_out.n_ov2 = cellfun(@(x,f) {f(x)},dis_rti,dis_info_in.n_ov2(req_pairs_i));
dis_info_out.cspace = dis_info_in.cspace(:,req_pairs_i);
dis_info_out.cspace_ex = dis_info_in.cspace_ex(:,req_pairs_i);

% dis_info_out.p_ov_ti=cellfun(@(x,i) {cast(ismemberloc(x(:,dis_version),i),'like',x)},dis_info_in.p_ov_ti(req_pairs_i),dis_rti);
% dis_info_out.p1_ov=cellfun(@(x) {x(:,:,dis_version)},dis_info_in.p1_ov(req_pairs_i));
% dis_info_out.p2_ov=cellfun(@(x) {x(:,:,dis_version)},dis_info_in.p2_ov(req_pairs_i));
dis_info_out.p1_ov_lim=cellfun(@(x,f) {f(x,:)},dis_rti,dis_info_in.p1_ov_lim(req_pairs_i));
dis_info_out.p2_ov_lim=cellfun(@(x,f) {f(x,:)},dis_rti,dis_info_in.p2_ov_lim(req_pairs_i));


dis_info_out.pics = full_pics(reduced_parts_idx);
dis_info_out.n_pics = numel(dis_info_out.pics);
dis_info_out.pairs = nchoosek(1:dis_info_out.n_pics,2);
dis_info_out.n_pairs = size(dis_info_out.pairs,1);
dis_info_out.non_ov_val = dis_info_in.non_ov_val(reduced_parts_idx,reduced_parts_idx);
dis_info_out.rot = dis_info_in.rot;
end
function loc=ismemberloc(a,b)
[~,loc]=ismember(a,b);
end
function i=sorted_ind(A,dim,dir)
[~,i]=sort(A,dim,dir);
end