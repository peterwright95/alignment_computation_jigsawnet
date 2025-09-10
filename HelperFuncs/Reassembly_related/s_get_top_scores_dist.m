function [ind]=s_get_top_scores_dist(scores_i,n_top_scores,iou_thresh,ov_lim_main,ov_lim_addt,prevind,tforms,close_thresh)
ind = zeros(n_top_scores,1);
i=1;
j=1;
% if ~exist('close_thresh','var') && numel(scores_i) ~= size(ov_lim_main,1)&& numel(scores_i) ~= size(ov_lim_addt,1)
%     error('number of overlap masks (main|addt : %d|%d) not matching number of scores (%d)',size(ov_lim_main,1),size(ov_lim_addt,1),numel(scores_i));
% end
while(i<=n_top_scores)
    if j>numel(scores_i)
        ind=ind(1:i-1);
        return;
    end
    basei=[tocolumn(ind(1:i-1));tocolumn(prevind)];
    if exist('close_thresh','var')
        is_j_valid = is_new_ind_valid_by_tformdist(scores_i(j),basei,tforms,close_thresh);
    else
        is_j_valid = is_new_ind_valid_by_ov_lim(scores_i(j),basei,ov_lim_main,iou_thresh) &&...
                     is_new_ind_valid_by_ov_lim(scores_i(j),basei,ov_lim_addt,iou_thresh);
    end
    if  is_j_valid
        ind(i)=scores_i(j);
        i = i+1;
    end
    j=j+1;
end
end
function [isvalid]=is_new_ind_valid_by_ov_lim(candi,basei,ov_lim,thresh)
    cand_ov_lim = ov_lim(candi,:);
    base_ov_lim = ov_lim(basei,:);
%     D=s_ov_lim_iou_distfun(cand_ov_lim,base_ov_lim);
    D=c_iou_intervals(base_ov_lim',cand_ov_lim');
    isvalid=~any(D>thresh);    
end
function [isvalid]=is_new_ind_valid_by_ov(candi,basei,ov,thresh)
    cand_ov = ov(candi,:);
    base_ov = ov(basei,:);
    ov_intersect = sum(bsxfun(@and,base_ov,cand_ov),2);
    ov_union = sum(bsxfun(@or,base_ov,cand_ov),2);
    isvalid=~any(ov_intersect./ov_union>thresh);    
end
function [isvalid]=is_new_ind_valid_by_tformdist(candi,basei,tforms,thresh)
    cand_tform = tforms(candi,:);
    base_tforms = tforms(basei,:);
    isvalid = ~any(s_tform_dist(base_tforms,cand_tform)<thresh);
end