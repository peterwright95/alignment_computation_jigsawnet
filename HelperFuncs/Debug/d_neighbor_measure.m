function [n,correct,total,avgdistance,avg_iou]=d_neighbor_measure(pics,absT,randrot,iou_thresh)
if ~exist('iou_thresh','var')
    iou_thresh = 0.65;
end
n=0;correct=0;total=0;avgdistance=0;avg_iou=0;
n_pics=numel(pics);
% load_options.is_use_rand_frag_rot = true;
% load_options.rand_angles = randrot;
load_options.reduce_poly_length_to = 1.0;
[f,~] = s_loadpics(pics, load_options);
pics_info = s_pics_info(pics);
parts=cellfun(@(x) str2double(x.part),pics_info);
pairs=nchoosek(1:n_pics,2);
n_pairs=size(pairs,1);
series=join(unique(cellfun(@(x) {x.series},pics_info)),'+');

% if numel(series)>1
%     warning('Could not calculate neighbor measure because there are %d groups of fragments',numel(series));
%     return;
% end
% [ ~,rel_gtT] = s_get_ground_truth(pics_info{1}.series,pairs);
% rel_gtT = zeros(n_pairs,3);
% found_gt = false(n_pairs,1);
% for i=1:n_pairs
%     p1series=pics_info{pairs(i,1)}.series;
%     p2series=pics_info{pairs(i,2)}.series;
%     p1part=pics_info{pairs(i,1)}.part;
%     p2part=pics_info{pairs(i,1)}.part;
%     if strcmp(p1series,p2series)
%         [ ~,rel_gtT_res] = s_get_ground_truth( p1series,[str2double(p1part),str2double(p2part)]);
%         if ~isempty(rel_gtT_res)
%             rel_gtT(i,:) = rel_gtT_res;
%             found_gt(i) = true;
%         end
%     end
% end
[ ~,rel_gtT] = s_get_ground_truth(series,parts(pairs));
if isempty(rel_gtT)
    fprintf('Could not calculate neighbor meausre due to missing ground truth\n');
    return;
end
tform_p1 = absT(pairs(:,1),:);
tform_p2 = absT(pairs(:,2),:);
rel_T = s_compose_tforms(s_invtform(tform_p1),tform_p2);

% spoly = arrayfun(@(i) {s_shift_poly(f{i}.poly,abs_gtT(i,:))},1:n_pics);
% spoly_ex = arrayfun(@(i) {s_shift_poly(f{i}.poly_ex,abs_gtT(i,:))},1:n_pics);
sz = s_calc_globalsz(f,pairs);
gsz=max(sz);
gsz=[gsz,gsz];
masks=cellfun(@(x) {s_canvasSize(x.mask,gsz)},f);
masks_ex=cellfun(@(x) {s_canvasSize(x.mask_ex,gsz)},f);
total=0;
correct=0;
% avgdistance=0;
distances_sum=0;
iou_sum=0;
for pi=1:n_pairs
    p1=pairs(pi,1);
    p2=pairs(pi,2);
%     p1m=masks{p1};
    p2m=masks{p2};
    p1_exm=masks_ex{p1};
%     p2_exm=masks_ex{p2};
    
    sp2m_gt=transform_mask(p2m,rel_gtT(pi,:));
    if sum(sum(p1_exm&sp2m_gt))>0% || sum(sum(p1m&transform_mask(p2_exm,rel_gtT(pi,:))))>0
        total=total+1;
        cand_rel_tform = s_compose_tforms([0,0,-randrot(p1)],s_compose_tforms(rel_T(pi,:),[0,0,randrot(p2)]));
        sp2m_cand=transform_mask(p2m,cand_rel_tform);
        tintersection = sum(sum(sp2m_cand&sp2m_gt));
        tunion = sum(sum(sp2m_cand|sp2m_gt));
        
        spoly = s_shift_poly(f{p2}.poly,cand_rel_tform);
        d=polygons_dist(f{p1}.poly,f{p1}.poly_ex,spoly);
        distances_sum = distances_sum + d;
        
        
%         figure(1);imshow(double(masks{p1})/2+double(sp2m_gt)/3);title('gt');
%         figure(2);imshow(double(masks{p1})/2+double(sp2m_cand)/3);title('assemb');
%         figure(3);imshow(double(sp2m_cand&sp2m_gt)/2+double(sp2m_cand|sp2m_gt)/2);title(sprintf('iou: %0.4f',tintersection/tunion));

%         figure(4);imshow(d_Tforms2im(absT,s_loadpics(pics, struct('is_use_rand_frag_rot',true,'rand_angles',randrot))));title('result assemb');
%         figure(5);imshow(d_Tforms2im(s_get_ground_truth(pics_info{1}.series,pairs),f));title('gt assemb');
%         figure(6);plot(f{p1}.poly(1,:),f{p1}.poly(2,:),f{p1}.poly_ex(1,:),f{p1}.poly_ex(2,:),spoly(1,:),spoly(2,:));
        tiou=tintersection/tunion;
        iou_sum=iou_sum+tiou;
        if tiou > iou_thresh
            correct = correct + 1;
        else
            tmp=1;
        end
    end
end
n=correct/total;
avgdistance = distances_sum/total;
avg_iou=iou_sum/total;
% figure(4);imshow(d_Tforms2im(absT,s_loadpics(pics, struct('is_use_rand_frag_rot',true,'rand_angles',randrot))));title(sprintf('result assemb. neighbor measure: %0.4f, correct: %d, total: %d', n, correct, total));
end
function tm = transform_mask(mask,T)
tm=s_imshift(imrotate(mask,T(3),'crop'),T(1:2));
end
function d=polygons_dist(pol1,pol1_ex,pol2)
% TOL = 1.0e-12;
in=c_inpoly(pol2',pol1_ex');
if any(in)
    dists = p_poly_dist(pol2(1,in),pol2(2,in),pol1(1,:),pol1(2,:));
    d=max(mean(dists),0);
%     if any(dists<-TOL)
%         d=min(dists);
%     else
%         d=max(dists);
%     end
else
%     d=inf;r
    dists = p_poly_dist(pol2(1,:),pol2(2,:),pol1(1,:),pol1(2,:));
    d=min(dists);
end
end