function [ conf,match_info,allconf,allmatch_info ] = confidence( dis_info,options,fr_input )
% closeness_threshold = options.closeness_threshold;
% ov_cutoff = options.ov_cutoff;
n_pairs = dis_info.n_pairs;
if ~isfield(options,'selected_pairs')
    selected_pairs=1:n_pairs;
else
    selected_pairs=options.selected_pairs;
end
if ~isfield(options,'save_top_conf_frag_num')
    save_top_conf_frag_num = [];
else
    save_top_conf_frag_num = options.save_top_conf_frag_num;    
end
% N = dis_info.n_pics;
dis_info = fix_input(dis_info);
n_top_scores=10;
conf = -inf(n_pairs,1);
match_info=cell(n_pairs,1);
allconf = -inf(n_pairs,n_top_scores*n_top_scores);
allmatch_info=cell(n_pairs,1);
tAll=0;
for pi=torow(selected_pairs) %for each pair
    tPairMeas = tic;
    p=dis_info.pairs(pi,:);
    pimask=false(size(dis_info.pairs));pimask(pi,:)=true;
    p1_mask=dis_info.pairs==p(1);% & pimask; % pimask - to use only local neighbors
    p2_mask=dis_info.pairs==p(2);% & pimask; % pimask - to use only local neighbors

    [~,sorted_dis1_i]=sort(dis_info.dis1{pi});
    [~,sorted_dis2_i]=sort(dis_info.dis2{pi});
    valid_tforms_i = find(dis_info.valid_tform_mask{pi});
    if isempty(valid_tforms_i)
        tAll=tAll+toc(tPairMeas);
        continue;
    end
    vsorted_dis1_i = sorted_dis1_i(ismember(sorted_dis1_i,valid_tforms_i));
    vsorted_dis2_i = sorted_dis2_i(ismember(sorted_dis2_i,valid_tforms_i));
    p1_ov_ti = s_get_top_scores_dist(vsorted_dis1_i,n_top_scores,options.iou_thresh,dis_info.p1_ov_lim{pi},dis_info.p2_ov_lim{pi},[]);
    p2_ov_ti = s_get_top_scores_dist(vsorted_dis2_i,n_top_scores,options.iou_thresh,dis_info.p2_ov_lim{pi},dis_info.p1_ov_lim{pi},[]);
    uind = unique([p1_ov_ti;p2_ov_ti]);

    up1_ov_lim = dis_info.p1_ov_lim{pi}(uind,:);
    up2_ov_lim = dis_info.p2_ov_lim{pi}(uind,:);
%     [ucand_dis_p1,ucand_dis_p1_i2,ucand_dis_p1_ti,lim_iou1,lim_iou_pi_loc1] = find_2best_dis(uind,pi,up1_ov_lim,up2_ov_lim,p1_mask,dis_info,options);
%     [ucand_dis_p2,ucand_dis_p2_i1,ucand_dis_p2_ti,lim_iou2,lim_iou_pi_loc2] = find_2best_dis(uind,pi,up2_ov_lim,up1_ov_lim,p2_mask,dis_info,options);
    [ucand_dis_p1,ucand_dis_p1_i2,ucand_dis_p1_ti,closeT_dis_deriv_p1,n_closeT_p1] = find_2best_dis(uind,pi,up1_ov_lim,up2_ov_lim,p1_mask,dis_info,options);
    [ucand_dis_p2,ucand_dis_p2_i1,ucand_dis_p2_ti,closeT_dis_deriv_p2,n_closeT_p2] = find_2best_dis(uind,pi,up2_ov_lim,up1_ov_lim,p2_mask,dis_info,options);
    if ~isequal(n_closeT_p1,n_closeT_p2)
        warning('number of close neighbors changed, %d, %d',n_closeT_p1,n_closeT_p2);
    end
    [mconf,miou1,miou2,conf1,conf2,ind1,ind2,cdis_p1_2b_i2,cdis_p1_2b_ti,cdis_p2_2b_i1,cdis_p2_2b_ti]=...
        calc_mutual_conf(ucand_dis_p1,ucand_dis_p1_i2,ucand_dis_p1_ti,dis_info.dis1{pi}(uind),up1_ov_lim,...
                         ucand_dis_p2,ucand_dis_p2_i1,ucand_dis_p2_ti,dis_info.dis2{pi}(uind),up2_ov_lim,options.iou_thresh);

    [mconf_s,conf_i]=sort(mconf,'descend');

    n_all_conf = min(size(allconf,2),numel(mconf_s));
    
    if ~isempty(save_top_conf_frag_num) && any(p==save_top_conf_frag_num)
        isSaveAndClose = true;
        show_top_conf_tforms(dis_info,fr_input,pi,options,uind,ind1,ind2,conf_i,cdis_p1_2b_i2,cdis_p1_2b_ti,cdis_p2_2b_i1,cdis_p2_2b_ti,mconf,conf1,conf2,1,isSaveAndClose);
    end
    
% DEBUG: show_top_conf_tforms(dis_info,fr_input,pi,options,uind,ind1,ind2,conf_i,cdis_p1_2b_i2,cdis_p1_2b_ti,cdis_p2_2b_i1,cdis_p2_2b_ti,mconf,conf1,conf2,1)
% DEBUG: d_show_ov_limits_diff(fr_input,dis_info,uind,[],lim_iou2,lim_iou_pi_loc2,p2_mask,pi,[3,4],2959,78)
    allconf(pi,1:n_all_conf)=mconf_s(1:n_all_conf);
    conf(pi)=allconf(pi,1);
    conf_i=conf_i(1:n_all_conf);
    
    [match_info{pi},allmatch_info{pi}] = fill_match_info(uind,ind1,ind2,miou1,miou2,conf1,conf2,p1_ov_ti,p2_ov_ti,conf_i,conf(pi),...
        closeT_dis_deriv_p1,closeT_dis_deriv_p2,n_closeT_p1);
    
    if min(match_info{pi}.i1,match_info{pi}.i2)>n_top_scores
        warning('top confidence locations above n_top_scores %d: [%d,%d]',n_top_scores,conf_i1_loc_from_top(1),conf_i1_loc_from_top(2));
    end 
    
    tAll=tAll+toc(tPairMeas);
end
fprintf('Confidence Performance: Avg calc time for each pair: %0.3f sec, Total: %0.3f sec\n',tAll/numel(selected_pairs),tAll);
end
function show_top_conf_tforms(dis_info,fr_input,pi,options,uind,ind1,ind2,conf_i,...
                              cdis_p1_2b_i2,cdis_p1_2b_ti,cdis_p2_2b_i1,cdis_p2_2b_ti,mconf,conf1,conf2,top_i,isSaveAndClose)
    if ~exist('isSaveAndClose','var')
        isSaveAndClose = false;
    end
    p=dis_info.pairs(pi,:);
    figh=figure(1);
    bst_ti=uind(ind1(conf_i(top_i)));
    dis=dis_info.dis1{pi}(bst_ti);
    if isfield(options,'reltforms') && pi<=numel(options.reltforms)
        tform = options.reltforms{pi}(bst_ti,:);
    else
        tform = dis_info.tforms{pi}(bst_ti,:);
    end
    if ~isempty(tform)
        subplot(2,2,1);
        imshow(s_imcomb(fr_input{p(1)}.rgb,fr_input{p(2)}.rgb,fr_input{p(1)}.mask,fr_input{p(2)}.mask,tform));
        title(sprintf('p1, pair: (%d,%d), ti: %d, dis: %0.4f',p(1),p(2),bst_ti,dis));
    end
    
    bst_p=cdis_p1_2b_i2(ind1(conf_i(top_i)));
    bst_ti=cdis_p1_2b_ti(ind1(conf_i(top_i)));
    if bst_ti~=0
        [~,pairi]=ismember(sort([bst_p,p(1)]),sort(dis_info.pairs,2),'rows');
        bst_pair=dis_info.pairs(pairi,:);
        dis_ti=[dis_info.dis1{pairi}(bst_ti),dis_info.dis2{pairi}(bst_ti)];
        dis=dis_ti(bst_pair==p(1));
        if isfield(options,'reltforms') && pairi<=numel(options.reltforms)
            tform = options.reltforms{pairi}(bst_ti,:);
        else
            tform = dis_info.tforms{pairi}(bst_ti,:);
        end
        if ~isempty(tform)
            subplot(2,2,2);imshow(s_imcomb(fr_input{bst_pair(1)}.rgb,fr_input{bst_pair(2)}.rgb,fr_input{bst_pair(1)}.mask,fr_input{bst_pair(2)}.mask,tform));
            title(sprintf('p1, best pair: (%d,%d), best ti: %d, dis: %0.4f',bst_pair(1),bst_pair(2),bst_ti,dis));
        end
    end
    
    bst_ti=uind(ind2(conf_i(top_i)));
    dis=dis_info.dis2{pi}(bst_ti);
    if isfield(options,'reltforms') && pi<=numel(options.reltforms)
        tform = options.reltforms{pi}(bst_ti,:);
    else
        tform = dis_info.tforms{pi}(bst_ti,:);
    end
    if ~isempty(tform)
        subplot(2,2,3);imshow(s_imcomb(fr_input{p(1)}.rgb,fr_input{p(2)}.rgb,fr_input{p(1)}.mask,fr_input{p(2)}.mask,tform));
        title(sprintf('p2, pair: (%d,%d), ti: %d, dis: %0.4f',p(1),p(2),bst_ti,dis));
    end
    
    bst_p=cdis_p2_2b_i1(ind2(conf_i(top_i)));
    bst_ti=cdis_p2_2b_ti(ind2(conf_i(top_i)));
    if bst_ti~=0
        [~,pairi]=ismember(sort([bst_p,p(2)]),sort(dis_info.pairs,2),'rows');
        bst_pair=dis_info.pairs(pairi,:);
        dis_ti=[dis_info.dis1{pairi}(bst_ti),dis_info.dis2{pairi}(bst_ti)];
        dis=dis_ti(bst_pair==p(2));
        if isfield(options,'reltforms') && pairi<=numel(options.reltforms)
            tform = options.reltforms{pairi}(bst_ti,:);
        else
            tform = dis_info.tforms{pairi}(bst_ti,:);
        end
        if ~isempty(tform)
            subplot(2,2,4);imshow(s_imcomb(fr_input{bst_pair(1)}.rgb,fr_input{bst_pair(2)}.rgb,fr_input{bst_pair(1)}.mask,fr_input{bst_pair(2)}.mask,tform));
            title(sprintf('p2, best pair: (%d,%d), best ti: %d, dis: %0.4f',bst_pair(1),bst_pair(2),bst_ti,dis));
        end
    end
    suptitle(sprintf('Pair(%d,%d), Conf: %0.4f, Conf1: %0.4f, Conf2: %0.4f',p(1),p(2),mconf(conf_i(top_i)),conf1(ind1(conf_i(top_i))),conf2(ind2(conf_i(top_i)))));
    if isSaveAndClose
        set(gcf, 'Position', get(0, 'Screensize'));
        outfile = fullfile(pwd,'Pictures','Results','Assembly_progress','assemb_prog.jpg');
        save_fig(outfile);
        close(figh);
    end
end

function dis_info = fix_input(dis_info)
if ~iscell(dis_info.dis1)
    dis_info.dis1=num2cell(dis_info.dis1,1);
end
if ~iscell(dis_info.dis2)
    dis_info.dis2=num2cell(dis_info.dis2,1);
end
if ~iscell(dis_info.tforms)
    dis_info.tforms=num2cell(dis_info.tforms,[1,2]);
end
if ~isfield(dis_info,'valid_tform_mask')
    dis_info.valid_tform_mask=cellfun(@(x) {true(size(x))},dis_info.dis1);
end
dis_info.tforms=reshape(dis_info.tforms,1,numel(dis_info.tforms));
if size(dis_info.dis1,1)~=1
    dis_info.dis1=dis_info.dis1';
end
if size(dis_info.dis2,1)~=1
    dis_info.dis2=dis_info.dis2';
end
% if size(dis_info.p_ov_ti,1)~=1
%     dis_info.p_ov_ti=dis_info.p_ov_ti';
% end
% if size(dis_info.p1_ov,1)~=1
%     dis_info.p1_ov=dis_info.p1_ov';
% end
% if size(dis_info.p2_ov,1)~=1
%     dis_info.p2_ov=dis_info.p2_ov';
% end
if size(dis_info.p1_ov_lim,1)~=1
    dis_info.p1_ov_lim=dis_info.p1_ov_lim';
end
if size(dis_info.p2_ov_lim,1)~=1
    dis_info.p2_ov_lim=dis_info.p2_ov_lim';
end
if size(dis_info.valid_tform_mask,1)~=1
    dis_info.valid_tform_mask=dis_info.valid_tform_mask';
end
end
function [cand_dis,cand_dis_i,cand_dis_ti,ov_iou,lim_iou] = find_2best_dis_old(pi,ti_ov,ti_q_ov,ti_ov_lim,ti_q_ov_lim,p_locmask,dis_info,options)
    ov_iou=[];    
    ti_ov = permute(ti_ov,[3,2,1]);
%     ti_q_ov = permute(ti_q_ov,[3,2,1]);
    n_uind = size(ti_ov,3);
    
    p_locmask(cellfun(@isempty,dis_info.p_ov_ti),:)=false;
%     pov_1=dis_info.p1_ov(p_locmask(:,1));
%     pov_2=dis_info.p2_ov(p_locmask(:,2));
    pov_lim_1=dis_info.p1_ov_lim(p_locmask(:,1));
    pov_lim_2=dis_info.p2_ov_lim(p_locmask(:,2));
    n_p1=numel(find(p_locmask(:,1)));
    povti=[dis_info.p_ov_ti(p_locmask(:,1)),dis_info.p_ov_ti(p_locmask(:,2))];
%     nc_pov = cellfun(@(x) size(x,2),[pov_1,pov_2]);
%     if any(nc_pov ~= size(ti_ov,2))
%         error('numer of contour elements should be equal.\nfrag from other pairs: %s, curr frag: %d',num2str(nc_pov,'%d ,'),size(ti_ov,2));
%     end
%     tic;
%     ov_iou=cellfun(@(ov) {permute(tform_ov(ov,ti_ov),[1,3,2])},[pov_1,pov_2]);
%     time1=toc;
%     tic;
%     tic;
%     lim_iou2=cellfun(@(ovlim) {pdist2(ovlim,ti_ov_lim,@s_ov_lim_iou_distfun)},[pov_lim_1,pov_lim_2]);
%     toc;
%     tic;
    lim_iou=cellfun(@(ovlim) {c_iou_intervals(ovlim',ti_ov_lim')},[pov_lim_1,pov_lim_2]);
%     toc;
%     for idx=1:numel(lim_iou)
%         if max(abs(lim_iou2{idx}(:)-lim_iou{idx}(:)))>1e-10
%             error('iou is different');
%         end
%     end
%     fprintf('ov-lim time diff: %0.3f sec\n',time1-toc);
%     tf_ov_count=tform_ov(precalc_data.pov_mat,ti_ov);
%     tf_ov_count=mat2cell(tf_ov_count,precalc_data.ne_pov)';
    if find(p_locmask(pi,:))==1
%         pi_pov_q=dis_info.p2_ov{pi};
        pi_plim_q=dis_info.p2_ov_lim{pi};
    else
%         pi_pov_q=dis_info.p1_ov{pi};
        pi_plim_q=dis_info.p1_ov_lim{pi};
    end
    pairs_ov_i = [find(p_locmask(:,1));find(p_locmask(:,2))];
%     tf_q_ov_iou=permute(tform_ov(pi_pov_q,ti_q_ov),[1,3,2]);
%     tmptf_q_lim_iou=pdist2(pi_plim_q,ti_q_ov_lim,@s_ov_lim_iou_distfun);
    tf_q_lim_iou=c_iou_intervals(pi_plim_q',ti_q_ov_lim');
    tf_ov_count_pi_ind = find(pairs_ov_i==pi);
%     ov_iou{tf_ov_count_pi_ind}(ov_iou{tf_ov_count_pi_ind}>options.iou_thresh & tf_q_ov_iou>options.iou_thresh)=0;
    lim_iou{tf_ov_count_pi_ind}(lim_iou{tf_ov_count_pi_ind}>options.iou_thresh & tf_q_lim_iou>options.iou_thresh)=0;
    
    get_ov_ind=@(ti,ovc,i) torow(ti(ti~=0 & ovc(:,i)>options.ov_cutoff));
    if any(p_locmask(:,1))
        p1dis=dis_info.dis1(p_locmask(:,1));
        p1povti=povti(1:n_p1);
        p1tf_ov_count=lim_iou(1:n_p1);%ov_iou(1:n_p1);
        cand_dis1=arrayfun(@(i) {cellfun(@(pair_dis1,ti,ovc) {torow(pair_dis1(get_ov_ind(ti,ovc,i)))},p1dis,p1povti,p1tf_ov_count)},tocolumn(1:n_uind));
        cand_dis1_ti=arrayfun(@(i) {cellfun(@(ti,ovc) {get_ov_ind(ti,ovc,i)},p1povti,p1tf_ov_count)},tocolumn(1:n_uind));
        [worstdis1,worstdis1_ti]=cellfun(@(pair_dis1,ti) infmax(pair_dis1,ti(ti~=0)),p1dis,p1povti);
    else
        cand_dis1=cell(n_uind,1);
        cand_dis1_ti=cell(n_uind,1);
        worstdis1=[];
        worstdis1_ti=[];
    end
    if any(p_locmask(:,2))
        p2dis=dis_info.dis2(p_locmask(:,2));
        p2povti=povti(n_p1+1:end);
        p2tf_ov_count=lim_iou(n_p1+1:end);%ov_iou(n_p1+1:end);
        cand_dis2=arrayfun(@(i) {cellfun(@(pair_dis2,ti,ovc) {torow(pair_dis2(get_ov_ind(ti,ovc,i)))},p2dis,p2povti,p2tf_ov_count)},tocolumn(1:n_uind)) ;
        cand_dis2_ti=arrayfun(@(i) {cellfun(@(ti,ovc) {get_ov_ind(ti,ovc,i)},p2povti,p2tf_ov_count)},tocolumn(1:n_uind));
        [worstdis2,worstdis2_ti]=cellfun(@(pair_dis2,ti) infmax(pair_dis2,ti(ti~=0)),p2dis,p2povti);        
    else
        cand_dis2=cell(n_uind,1);
        cand_dis2_ti=cell(n_uind,1);
        worstdis2=[];
        worstdis2_ti=[];
    end
    cand_dis=cellfun(@(cdis1,cdis2) {[cdis1,cdis2]},cand_dis1,cand_dis2);
    
    ind=[dis_info.pairs(p_locmask(:,1),2);dis_info.pairs(p_locmask(:,2),1)]';
    cand_dis_i = cellfun(@(cdis) {cell2mat(cellfun(@(c,x) {repmat(c,1,numel(x))},num2cell(ind),cdis))},cand_dis);
    cand_dis = cellfun(@(cdis) {cell2mat(cdis)},cand_dis);
    cand_dis_ti=cellfun(@(cdis1ti,cdis2ti) {cell2mat([cdis1ti,cdis2ti])},cand_dis1_ti,cand_dis2_ti);
    worstdis = [worstdis1,worstdis2];
    worstdis_ti=[worstdis1_ti,worstdis2_ti];
    
    cand_dis_ti = cellfun(@(cdis,cdisti) {cdisti(~isinf(cdis))},cand_dis,cand_dis_ti);
    cand_dis_i = cellfun(@(cdis,cdisi) {cdisi(~isinf(cdis))},cand_dis,cand_dis_i);
    cand_dis = cellfun(@(cdis) {cdis(~isinf(cdis))},cand_dis);
    
%     empty_cdis = cellfun(@isempty,cand_dis);
%     if any(empty_cdis)
%         [minv,mini]=min(worstdis);
%         cand_dis(empty_cdis) = {minv};
%         cand_dis_i(empty_cdis) = {ind(mini)};
%         cand_dis_ti(empty_cdis) = {worstdis_ti(mini)};
%     end
end
function [cand_dis,cand_dis_i,cand_dis_ti,closeT_dis_deriv,n_closeT,lim_iou,lim_iou_pi_loc] = find_2best_dis(uind,pi,ti_ov_lim,ti_q_ov_lim,p_locmask,dis_info,options)
    pi_p_mask = p_locmask(pi,:);
    p_locmask(pi,:)=false;
    pi_ov_lim = cat(3,dis_info.p1_ov_lim{pi},dis_info.p2_ov_lim{pi});
    plim_iou = c_iou_intervals(pi_ov_lim(:,:,pi_p_mask)',ti_ov_lim',options.gamma_L);
    qlim_iou = c_iou_intervals(pi_ov_lim(:,:,~pi_p_mask)',ti_q_ov_lim',options.gamma_L);
    closeT_mask = plim_iou >= options.gamma_H & qlim_iou >= options.gamma_H;
    pi_pass_iou_mask = ~closeT_mask & plim_iou > options.gamma_L;
%     pi_pass_iou_mask = plim_iou > options.gamma_L;
%     pi_pass_iou_mask(sub2ind(size(pi_pass_iou_mask),uind,(1:numel(uind))'))=false;
    pidis=[dis_info.dis1{pi},dis_info.dis2{pi}];
    pdis=pidis(:,pi_p_mask);
    q=repmat(dis_info.pairs(pi,~pi_p_mask),numel(pdis),1);
    piti=(1:numel(pdis))';
    
    closeT_dis_deriv = arrayfun(@(i) dis2iou_derivative(plim_iou(:,i),qlim_iou(:,i),pdis,closeT_mask(:,i)),1:size(closeT_mask,2));
    n_closeT = sum(closeT_mask,1);
    
    all_pass_mask = pi_pass_iou_mask;
    all_dis = pdis;
    all_q = q;
    all_ti = piti;
%     lim_iou = {plim_iou};
    lim_iou = plim_iou;
    lim_iou_pi_loc = [1,size(plim_iou,1)];
    if any(p_locmask(:))
        other_pov_lim = [dis_info.p1_ov_lim(p_locmask(:,1))';dis_info.p2_ov_lim(p_locmask(:,2))'];
        other_lim_iou = c_iou_intervals(cell2mat(other_pov_lim)',ti_ov_lim',options.gamma_L);
        other_pass_iou_mask = other_lim_iou > options.gamma_L;
%         other_lim_iou=cellfun(@(ovlim) {c_iou_intervals(ovlim',ti_ov_lim',options.gamma_L)},other_pov_lim);
%         other_pass_iou_mask = cellfun(@(iou) {iou > options.gamma_L},other_lim_iou);
        other_dis=[dis_info.dis1(p_locmask(:,1))';dis_info.dis2(p_locmask(:,2))'];
        other_q=num2cell([dis_info.pairs(p_locmask(:,1),2);dis_info.pairs(p_locmask(:,2),1)]);
        other_q=cellfun(@(q,dis) {repmat(q,numel(dis),1)},other_q,other_dis);
        other_ti=cellfun(@(dis) {(1:numel(dis))'},other_dis);
        all_pass_mask = [all_pass_mask;other_pass_iou_mask];
%         all_pass_mask = [all_pass_mask;cell2mat(other_pass_iou_mask)];
        all_dis = [all_dis;cell2mat(other_dis)];
        all_q = [all_q;cell2mat(other_q)];
        all_ti = [all_ti;cell2mat(other_ti)];
        if nargout>5
            lim_iou = [lim_iou;other_lim_iou];
            n_other_poc_lim = cellfun(@(x) size(x,1),other_pov_lim)';
            st = cumsum([1,size(plim_iou,1),n_other_poc_lim]);
            lim_iou_pi_loc = [st(1:end-1)',(st(2:end)-1)'];
        end
    end
%     all_pass_maskcell=num2cell(all_pass_mask,1);
%     cand_dis_2 = cellfun(@(pass_mask) {all_dis(pass_mask)},all_pass_maskcell);
%     cand_dis_i_2 = cellfun(@(pass_mask) {all_q(pass_mask)},all_pass_maskcell);
%     cand_dis_ti_2 = cellfun(@(pass_mask) {all_ti(pass_mask)},all_pass_maskcell);
    [all_pass_mask_r,all_pass_mask_c]=find(all_pass_mask);
    cand_dis = accumarray(all_pass_mask_c,all_dis(all_pass_mask_r),[size(all_pass_mask,2),1],@(x){x},[])';
    cand_dis_i = accumarray(all_pass_mask_c,all_q(all_pass_mask_r),[size(all_pass_mask,2),1],@(x){x},[])';
    cand_dis_ti = accumarray(all_pass_mask_c,all_ti(all_pass_mask_r),[size(all_pass_mask,2),1],@(x){x},[])';
end
function d=dis2iou_derivative(piou,qiou,fdis,mask)
iou=mean([piou(mask),qiou(mask)],2);
dis=fdis(mask);
d=sum(dis.*(iou./sum(iou)));
end
function count=tform_ov(ov,compared_ov_row)
    ov_intersect=sum(bsxfun(@and,ov,compared_ov_row),2);
    ov_union=sum(bsxfun(@or,ov,compared_ov_row),2);
    count = ov_intersect ./ ov_union;
end
function [mconf,D1,D2,conf1,conf2,ind1,ind2,cdis_p1_2b_i2,cdis_p1_2b_ti,cdis_p2_2b_i1,cdis_p2_2b_ti]=...
        calc_mutual_conf(ucdis_p1,ucdis_p1_i2,ucdis_p1_ti,dis1,up1_ov_lim,ucdis_p2,ucdis_p2_i1,ucdis_p2_ti,dis2,up2_ov_lim,iou_thresh)
    if ~iscell(ucdis_p1)
        error('ucdis_p1 not cell. content=%s',mat2str(ucdis_p1));
    end
    [ucdis_p1_s,cdis_p1_si]=cellfun(@(x) sort(x,'ascend'),ucdis_p1,'UniformOutput',false);
    empty_cdis_p1_mask = cellfun(@isempty,ucdis_p1);
    cdis_p1_2b=zeros(size(dis1));cdis_p1_2b_i2=zeros(size(dis1));cdis_p1_2b_ti=zeros(size(dis1));
    cdis_p1_2b(~empty_cdis_p1_mask)=cellfun(@(x) x(1),ucdis_p1_s(~empty_cdis_p1_mask));
    cdis_p1_2b_i2(~empty_cdis_p1_mask)=cellfun(@(i,x) i(x(1)),ucdis_p1_i2(~empty_cdis_p1_mask),cdis_p1_si(~empty_cdis_p1_mask));
    cdis_p1_2b_ti(~empty_cdis_p1_mask)=cellfun(@(i,x) i(x(1)),ucdis_p1_ti(~empty_cdis_p1_mask),cdis_p1_si(~empty_cdis_p1_mask));
    
    [ucdis_p2_s,cdis_p2_si]=cellfun(@(x) sort(x,'ascend'),ucdis_p2,'UniformOutput',false);
    empty_cdis_p2_mask = cellfun(@isempty,ucdis_p2);
    cdis_p2_2b=zeros(size(dis2));cdis_p2_2b_i1=zeros(size(dis2));cdis_p2_2b_ti=zeros(size(dis2));
    cdis_p2_2b(~empty_cdis_p2_mask)=cellfun(@(x) x(1),ucdis_p2_s(~empty_cdis_p2_mask));
    cdis_p2_2b_i1(~empty_cdis_p2_mask)=cellfun(@(i,x) i(x(1)),ucdis_p2_i1(~empty_cdis_p2_mask),cdis_p2_si(~empty_cdis_p2_mask));
    cdis_p2_2b_ti(~empty_cdis_p2_mask)=cellfun(@(i,x) i(x(1)),ucdis_p2_ti(~empty_cdis_p2_mask),cdis_p2_si(~empty_cdis_p2_mask));
    
    conf1=zeros(size(dis1));conf2=zeros(size(dis2));
    conf1(~empty_cdis_p1_mask)=(1-dis1(~empty_cdis_p1_mask)./(cdis_p1_2b(~empty_cdis_p1_mask)+eps));
    conf2(~empty_cdis_p2_mask)=(1-dis2(~empty_cdis_p2_mask)./(cdis_p2_2b(~empty_cdis_p2_mask)+eps));
        
    [ind1,ind2]=meshgrid(1:numel(conf1),1:numel(conf2));
    ind1=ind1(:);ind2=ind2(:);
    mconf = (conf1(ind1)+conf2(ind2))./2;
%     mconf = max(conf1(ind1),conf2(ind2));
    
    [isvalid,D1,D2]=is_new_ind_valid_by_ov_lim(ind1,ind2,up1_ov_lim,up2_ov_lim,iou_thresh);
    mconf(~isvalid) = -inf;
end
function [isvalid,iou1,iou2]=is_tf_pair_valid_by_ov(ind1,ind2,up1_ov,up2_ov,iou_thresh)
    p2_ov_pi_i1=up2_ov(ind1,:);
    p2_ov_pi_i2=up2_ov(ind2,:);
    inters_count = sum(p2_ov_pi_i1 & p2_ov_pi_i2,2);
    union_count = sum(p2_ov_pi_i1 | p2_ov_pi_i2,2);
    iou2 = inters_count./union_count;
    p1_ov_pi_i1=up1_ov(ind1,:);
    p1_ov_pi_i2=up1_ov(ind2,:);
    inters_count = sum(p1_ov_pi_i1 & p1_ov_pi_i2,2);
    union_count = sum(p1_ov_pi_i1 | p1_ov_pi_i2,2);
    iou1 = inters_count./union_count;
    isvalid = iou2 >= iou_thresh & iou1 >= iou_thresh;
end
function [isvalid,iou1,iou2]=is_new_ind_valid_by_ov_lim(ind1,ind2,up1_ov_lim,up2_ov_lim,iou_thresh)
    D1=c_iou_intervals(up1_ov_lim',up1_ov_lim');
    D2=c_iou_intervals(up2_ov_lim',up2_ov_lim');
    ind = sub2ind([size(up1_ov_lim,1),size(up1_ov_lim,1)],ind1,ind2);
    iou1=D1(ind);
    iou2=D2(ind);
    isvalid=iou1 >= iou_thresh & iou2 >= iou_thresh;
end
function [match_info,allmatch_info] = fill_match_info(uind,ind1,ind2,miou1,miou2,conf1,conf2,p1_ov_ti,p2_ov_ti,conf_i,piconf,...
    p1deriv,p2deriv,ncand)
    allmatch_info.tform_i1=uind(ind1(conf_i));
    allmatch_info.tform_i2=uind(ind2(conf_i));
    allmatch_info.iou1=miou1(conf_i);
    allmatch_info.iou2=miou2(conf_i);
    allmatch_info.conf1=conf1(ind1(conf_i));
    allmatch_info.conf2=conf2(ind2(conf_i));
    allmatch_info.i1 = func_output(2,@ismember,uind(ind1(conf_i)),p1_ov_ti);
    allmatch_info.i2 = func_output(2,@ismember,uind(ind2(conf_i)),p2_ov_ti);
    allmatch_info.deriv1=p1deriv(ind1(conf_i));
    allmatch_info.deriv2=p2deriv(ind2(conf_i));
    allmatch_info.ncand=ncand(ind1(conf_i));
        
    match_info.tform_i1=allmatch_info.tform_i1(1);
    match_info.tform_i2=allmatch_info.tform_i2(1);
    match_info.iou1=allmatch_info.iou1(1);
    match_info.iou2=allmatch_info.iou2(1);
    match_info.conf1=allmatch_info.conf1(1);
    match_info.conf2=allmatch_info.conf2(1);
    match_info.i1 = allmatch_info.i1(1);
    match_info.i2 = allmatch_info.i2(1);
    match_info.deriv1=allmatch_info.deriv1(1);
    match_info.deriv2=allmatch_info.deriv2(1);
    match_info.ncand=allmatch_info.ncand(1);
    match_info.debug_str = {sprintf("conf1:%0.4f, iou1:%0.4f\ni1:%d,conf:%0.4f",...
        match_info.conf1,match_info.iou1,match_info.i1,piconf),...
                                     sprintf("conf2:%0.4f, iou2:%0.4f\ni2:%d,conf:%0.4f",...
        match_info.conf2,match_info.iou2,match_info.i2,piconf)};
end