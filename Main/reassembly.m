function [ fullrgb, T, R,order,totalconf,Adj,fullmask ] = reassembly(input, fr_info, dis_info, options)

POLYS_CLOSENESS_THRESHOLD = 30;
N = numel(input);
Q = 1:N;
P = {[]};
T = zeros(N,2);
R = zeros(N,1);
fullrgb = {[]};
fullmask = {[]};
full_i = 1;
Adj = {[]};%zeros(N);
im = cellfun(@(x) {x.rgb},input);
masks = cellfun(@(x) {x.mask},input);
masks_ex = cellfun(@(x) {x.mask_ex},input);
rel_tforms = dis_info.tforms;
abs_tforms = cell(size(rel_tforms));%zeros(size(rel_tforms));
is_valid_tforms = cellfun(@(x) {true(size(x,1),1)},rel_tforms);%true(size(rel_tforms,1),size(dis_info.pairs,1));
calc_pairs_tforms_mask = false(size(dis_info.pairs,1),1);
totalconf=[];
final_polys = cell(N,1);
% mft_cs = cell(N,N,numel(dis_info.rot)); %multi fragment tform c-space
P_info=[];
pics_info=s_pics_info(dis_info.pics);
Adj{full_i}=zeros(N);
totalconf(full_i) = 0;

while ~isempty(Q)
    fprintf('Reassemb in progress... Fragments placed: %d, Fragments remaining: %d\n',numel(P{full_i}),numel(Q));
    if (~isempty(P{full_i}))
        lia = ismember(dis_info.pairs,P{full_i});
        [abs_tforms,calc_pairs_tforms_mask] = get_abs_tforms(rel_tforms,abs_tforms,lia,T,R,dis_info.pairs,calc_pairs_tforms_mask);
        [is_valid_tforms] = disable_tforms_cspace(abs_tforms,is_valid_tforms,[T,R],lia,dis_info.pairs,dis_info.cspace,Q);
%         [multiq_tforms,multiqtf_p]=calc_multi_frag_tform(lia,dis_info.pairs,P{full_i},Q,dis_info.rot,options.sampling_res,[T,R], fullmask{full_i},masks,masks_ex,dis_info.cspace,dis_info.cspace_ex,is_valid_tforms);
        [multiq_tforms,multiqtf_p]=recalc_multifr_tforms(lia,dis_info.pairs,P{full_i},Q,[T,R],dis_info.cspace_ex,is_valid_tforms,abs_tforms);
        [conf2,confP,confq,bst_tforms,Q2Pmap]=reassemb_confidence(P_info,abs_tforms,is_valid_tforms,dis_info,P{full_i},Q,input,T,R,fr_info,multiq_tforms,multiqtf_p,options);
        [top_i,qi,maxconf,qtform,bb_mask]=nextp(conf2,confP,confq,bst_tforms);

%         d_show_reassmb_progress( conf2,P{full_i},Q,T,R,bst_tforms,input,top_i,qi);
        [~,pairs_with_q_i]=find(dis_info.pairs'==Q(qi));
        pairs_with_q_mask = ismember(dis_info.pairs(pairs_with_q_i,:),Q);
        pairs_with_q_not_in_P_i = pairs_with_q_i(sum(pairs_with_q_mask,2)==2);
        [sorted_Q_pair_conf,sorted_Q_pair_conf_i] = sort(dis_info.conf(pairs_with_q_not_in_P_i),'descend');
        
        if (numel(sorted_Q_pair_conf) < 2 || ... %dont open new image if only 2 fragments in Q (only 1 or 0 pairs with 1)
            maxconf >= sorted_Q_pair_conf(2) ||...
           full_i == options.max_output_images)
            q=Q(qi);
            fprintf('fragment %d was chosen as next piece with conf: %0.3f. Matched with fragments: %s\n',...
                q,maxconf,mat2str(Q2Pmap{qi,top_i}));
%             p=Q2Pmap{qi,top_i}(1);
            totalconf(full_i)=totalconf(full_i)+maxconf;
            Adj{full_i}(Q2Pmap{qi,top_i}, repmat(q,numel(Q2Pmap{qi,top_i}),1)) = maxconf;
            Adj{full_i}(repmat(q,numel(Q2Pmap{qi,top_i}),1),Q2Pmap{qi,top_i}) = maxconf;
        else
            pair_pq_i = pairs_with_q_not_in_P_i(sorted_Q_pair_conf_i(1));
            p=pairs(pair_pq_i,1);
            q=pairs(pair_pq_i,2);
            fprintf('starting a new image with pair %d,%d\n',p,q);
            full_i = full_i + 1;
            Adj{full_i}=zeros(N);
            totalconf(full_i) = 0;
            
            totalconf(full_i)=totalconf(full_i)+sorted_Q_pair_conf(1);

            qtform = s_tform_mean([dis_info.tforms{pair_pq_i}(dis_info.best_match_info{pair_pq_i}.tform_i1,:);...
                                   dis_info.tforms{pair_pq_i}(dis_info.best_match_info{pair_pq_i}.tform_i2,:)]);
            

            Q(Q==p) = [];
            P{full_i}=p;
            Adj{full_i} = zeros(N);
            Adj{full_i}(p, q) = sorted_Q_pair_conf(1);
            Adj{full_i}(q, p) = sorted_Q_pair_conf(1);
            fullrgb{full_i} = [];
        end
    else
        [sorted_conf_val,sorted_conf_i]=sort(dis_info.conf,'descend');
        starting_pair_i=sorted_conf_i(options.start_with_pair);
        
        totalconf(full_i)=totalconf(full_i)+sorted_conf_val(options.start_with_pair);

        qtform = s_tform_mean([dis_info.tforms{starting_pair_i}(dis_info.best_match_info{starting_pair_i}.tform_i1,:);...
                               dis_info.tforms{starting_pair_i}(dis_info.best_match_info{starting_pair_i}.tform_i2,:)]); 
        p=dis_info.pairs(starting_pair_i,1);
        q=dis_info.pairs(starting_pair_i,2);
        Q(Q==p) = [];
        P{full_i}=[P{full_i},p];
        
        Adj{full_i}(p, q) = sorted_conf_val(options.start_with_pair);
        Adj{full_i}(q, p) = sorted_conf_val(options.start_with_pair);

    end
    fprintf('fragment %d absolute selected sampled tform (before refine): %s\n',q,mat2str(qtform));
    out_candT = qtform;
%     [ out_candT,fval ] = s_refine_w_measure( P_info,fr_info,input,P{full_i},q,[T,R],qtform,options);
    [out_candT,out_fval] = refine_q_tform( P_info,fr_info,input,P{full_i},q,[T,R],qtform,dis_info.pairs,dis_info.cspace,options );
    fprintf('fragment %d absolute selected sampled tform (after refine): %s\n',q,mat2str(out_candT));
%     fprintf('refine on frag %d changed the tform by %0.3f (dist)\n',q,s_tform_dist([t,r],out_candT));
    Q(Q==q) = [];
    P{full_i}=[P{full_i},q];
    T(q,:) = out_candT(1:2);
    R(q) = out_candT(3);
    mq = im{q};
    
    if (isempty(fullrgb{full_i}))
        fullrgb{full_i} = im{p};
        fullmask{full_i} = masks{p};
        final_polys{p} = s_mask2poly(masks{p});
    end
    [uP_info,Debug_masks]=find_P_mask_order(input,fr_info,P_info,P{full_i},q,[T,R],options);
    if isempty(uP_info) % assume a new fragment was selected to be placed inside the P polygon. restore previous info.
        uP_info = P_info;
    end
    %     figure(1);subplot(1,2,1);imshow(Debug_masks{1});subplot(1,2,2);imshow(Debug_masks{2});figure(2);imshow(convert_distance_color(Debug_masks{3}));colormap(jet);colorbar;
    %     figure(1);subplot(1,2,1);imshow(Debug_masks{4});subplot(1,2,2);plot([uP_info.poly(1,:),uP_info.poly(1,1)],[uP_info.poly(2,:),uP_info.poly(2,1)]);
    bq = masks{q};
    [noref_res_rgb,noref_res_mask] = s_imcomb(fullrgb{full_i},mq,fullmask{full_i}, bq,qtform);
    [wref_res_rgb,wref_res_mask,~,bq_tformed] = s_imcomb(fullrgb{full_i},mq,fullmask{full_i}, bq,out_candT);
    %     figure(1);subplot(1,2,1);imshow(s_imtrim(wref_res_rgb,wref_res_mask));title(['after - ' mat2str(out_candT)]);subplot(1,2,2);imshow(s_imtrim(noref_res_rgb,noref_res_mask));title(['before' mat2str(qtform)]);
    if options.is_output_assemb_progress
        progress_wref = s_imtrim(wref_res_rgb,wref_res_mask);
        progress_woref = s_imtrim(noref_res_rgb,noref_res_mask);
        progsz = max(size(progress_woref),size(progress_wref));
        progress_woref = insertText(s_canvasSize(progress_woref,progsz),[2,2],['before refinement on fr #' num2str(q)],'FontSize',20);
        progress_wref = insertText(s_canvasSize(progress_wref,progsz),[2,2],['after refinement on fr #' num2str(q)],'FontSize',20);
        prog_out_filename = fullfile(options.progress_out_folder,[pics_info{q}.series '_assemb_progress.jpg']);
        uimwrite([progress_wref,progress_woref],prog_out_filename);
        
        if ~isempty(Debug_masks{1}) && ~isempty(Debug_masks{2})
            dsz=max(size(Debug_masks{1}),size(Debug_masks{2}));
            prog_debug_out_filename = fullfile(options.progress_out_folder,[pics_info{q}.series '_assemb_progress_boundary.jpg']);
            uimwrite([s_canvasSize(Debug_masks{1},dsz),s_canvasSize(Debug_masks{2},dsz)],prog_debug_out_filename);
        end
    end
    P_info = uP_info;
    
    fullrgb{full_i}=wref_res_rgb;
    fullmask{full_i}=wref_res_mask;
    options.curr_assemb = fullrgb{full_i};
    options.curr_assemb_mask = fullmask{full_i};
    final_polys{q} = s_mask2poly(bq_tformed);
    not_connected2q = intersect(find(Adj{full_i}(q,:)==0),P{full_i}(P{full_i}~=q));
    for cand = not_connected2q
        dist = abs(s_polys_dist(final_polys{cand},final_polys{q}));
        if dist<POLYS_CLOSENESS_THRESHOLD 
            Adj{full_i}(cand,q) = -1;
            Adj{full_i}(q,cand) = -1;
        end
    end
end
for i=1:full_i
    [fullrgb{i},~,fullmask{i}] = s_imtrim(fullrgb{i},fullmask{i});
end
order=P;
fprintf('\n!!!Re-assembly calculation finished!!!\n');
end
function [P_info,Debug_masks]=find_P_mask_order(fr_input_info,fr_info,P_info,Pq,q,TR,options)
translate_points=@(points_rc,tr) round(bsxfun(@plus,points_rc,tr));
% shift_mask2rc=@(fi,mask,gsz) translate_points(s_rotate_mask_with_ind(find(s_canvasSize(mask,gsz)),gsz,TR(fi,3)),TR(fi,1:2));
shift_mask2rc=@(fi,mask,gsz) translate_points(s_rotate_ind_loose( mask_indices(mask),size(mask),gsz,TR(fi,3) ), TR(fi,1:2));
p=Pq(Pq~=q);
globalsz=s_calclimits(fr_input_info,fr_info,TR,p,q);
Debug_masks=cell(1,4);
if isempty(P_info)
    if numel(Pq)==2
        sp=shift_mask2rc(p,fr_info{p}.c_mask,globalsz);
        sp_ex=shift_mask2rc(p,fr_info{p}.ex_only_mask,globalsz);
        pMi=fr_input_info{p}.mask_Mi;
        pSi=fr_input_info{p}.mask_Si;
        p_ex_only_vals=fr_info{p}.ex_only_vals;
        p_ex_only_vals_rgb=fr_info{p}.ex_only_vals_rgb;
        p_c_vals=fr_info{p}.c_vals;
        p_c_vals_rgb=fr_info{p}.c_vals_rgb;
        p_g_info=s_update_g_info_rot(fr_info{p}.g_info,TR(p,3));
    else
        error('initialization should happen when numel(P)==2. instead was %d',numel(Pq));
    end
else
    sp=P_info.P_rc; %adjust to globalsz
    sp_ex=P_info.P_ex_rc; %adjust to globalsz
    if any(globalsz>P_info.globalsz)
       offset=round((globalsz-P_info.globalsz)./2);
       sp=bsxfun(@plus,sp,offset);
       sp_ex=bsxfun(@plus,sp_ex,offset);
    else
        globalsz=P_info.globalsz;
    end
    pMi=P_info.PMi;
    pSi=P_info.PSi;
    p_ex_only_vals=P_info.P_ex_vals;
    p_ex_only_vals_rgb=P_info.P_ex_vals_rgb;
%     p_ex_gdir=P_info.P_ex_gdir;
    p_c_vals=P_info.P_vals;
    p_c_vals_rgb=P_info.P_vals_rgb;
%     p_c_gdir=P_info.P_gdir;
%     p_c_gw8s=P_info.P_g_w8s;
    p_g_info=P_info.g_info;
end
P_info=[];
sq=shift_mask2rc(q,fr_info{q}.c_mask,globalsz);
sq_ex=shift_mask2rc(q,fr_info{q}.ex_only_mask,globalsz);
ov_on_p=ismember(sp,sq_ex,'rows');
ov_on_q=ismember(sq,sp_ex,'rows');
qMi=fr_input_info{q}.mask_Mi;
qSi=fr_input_info{q}.mask_Si;
ov_lim_on_p = s_ov2limits(ov_on_p,pMi);
ov_lim_on_q = s_ov2limits(ov_on_q,qMi);
if any(isnan(ov_lim_on_p))||any(isnan(ov_lim_on_q))
    warning('No overlap detected for final transformation');
    return;
end
[~,ploc]=ismember(ov_lim_on_p,pMi);
[~,qloc]=ismember(ov_lim_on_q,qMi);
ov_lim_on_p_rc = sp(ploc,:);
ov_lim_on_q_rc = sq(qloc,:);

%assume p1 is the start point
plimSi=round(ov_lim_on_p*numel(pMi));
qlimSi=round(ov_lim_on_q*numel(qMi));
if ov_lim_on_q(1)<ov_lim_on_q(2)
    q1_q2=[sq(qSi(qlimSi(1)-1:-1:1),:);sq(qSi(end:-1:qlimSi(2)+1),:)];
else
    q1_q2=sq(qSi(qlimSi(1)-1:-1:qlimSi(2)+1),:);
end
if ov_lim_on_p(1)<ov_lim_on_p(2)
    p2_p1=[sp(pSi(plimSi(2)+1:end),:);sp(pSi(1:plimSi(1)-1),:)];
else
    p2_p1=sp(pSi(plimSi(2)+1:plimSi(1)-1),:);
end
line_p1q1_x=[ov_lim_on_p_rc(1,2),ov_lim_on_q_rc(1,2)];
line_p1q1_y=[ov_lim_on_p_rc(1,1),ov_lim_on_q_rc(1,1)];
line_p2q2_x=[ov_lim_on_p_rc(2,2),ov_lim_on_q_rc(2,2)];
line_p2q2_y=[ov_lim_on_p_rc(2,1),ov_lim_on_q_rc(2,1)];
line_p1q1_norm=sqrt((line_p1q1_x(1)-line_p1q1_x(2)).^2+(line_p1q1_y(1)-line_p1q1_y(2)).^2);
line_p2q2_norm=sqrt((line_p2q2_x(1)-line_p2q2_x(2)).^2+(line_p2q2_y(1)-line_p2q2_y(2)).^2);
norm_p1q1_p2q2=line_p1q1_norm+line_p2q2_norm;
is_intersect_p1q1_p2q2=~isempty(polyxpoly(line_p1q1_x,line_p1q1_y,line_p2q2_x,line_p2q2_y));
line_p1q2_x=[ov_lim_on_p_rc(1,2),ov_lim_on_q_rc(2,2)];
line_p1q2_y=[ov_lim_on_p_rc(1,1),ov_lim_on_q_rc(2,1)];
line_p2q1_x=[ov_lim_on_p_rc(2,2),ov_lim_on_q_rc(1,2)];
line_p2q1_y=[ov_lim_on_p_rc(2,1),ov_lim_on_q_rc(1,1)];
line_p1q2_norm=sqrt((line_p1q2_x(1)-line_p1q2_x(2)).^2+(line_p1q2_y(1)-line_p1q2_y(2)).^2);
line_p2q1_norm=sqrt((line_p2q1_x(1)-line_p2q1_x(2)).^2+(line_p2q1_y(1)-line_p2q1_y(2)).^2);
norm_p1q2_p2q1=line_p1q2_norm+line_p2q1_norm;
is_intersect_p1q2_p2q1=~isempty(polyxpoly(line_p1q2_x,line_p1q2_y,line_p2q1_x,line_p2q1_y));
if is_intersect_p1q2_p2q1 || (~is_intersect_p1q2_p2q1 && ~is_intersect_p1q1_p2q2 &&...
        norm_p1q1_p2q2<norm_p1q2_p2q1)
    [p1_q1(:,2),p1_q1(:,1)] = line2mask(ov_lim_on_p_rc(1,2),ov_lim_on_p_rc(1,1),ov_lim_on_q_rc(1,2),ov_lim_on_q_rc(1,1));
    [p2_q2(:,2),p2_q2(:,1)] = line2mask(ov_lim_on_p_rc(2,2),ov_lim_on_p_rc(2,1),ov_lim_on_q_rc(2,2),ov_lim_on_q_rc(2,1));
    orderedP_rc=[p1_q1;q1_q2;flip(p2_q2,1);p2_p1];            
else
    [p1_q2(:,2),p1_q2(:,1)] = line2mask(ov_lim_on_p_rc(1,2),ov_lim_on_p_rc(1,1),ov_lim_on_q_rc(2,2),ov_lim_on_q_rc(2,1));
    [p2_q1(:,2),p2_q1(:,1)] = line2mask(ov_lim_on_p_rc(2,2),ov_lim_on_p_rc(2,1),ov_lim_on_q_rc(1,2),ov_lim_on_q_rc(1,1));
    orderedP_rc=[p1_q2;flip(q1_q2,1);flip(p2_q1,1);p2_p1];
end
[~,~,orderedP_rc_ui]=unique(orderedP_rc,'rows','stable');
n_orderedP = numel(orderedP_rc_ui);
orderedP_rc_u_loc=accumarray(orderedP_rc_ui,(1:n_orderedP)',[],@(x){x});
n_orderedP_rc_ui=cellfun(@numel,orderedP_rc_u_loc)>1;
% P_no_loop_mask=cat(2,true(n_orderedP,1),...
%                cell2mat(cellfun(@(x) {[true(min(x),1);false(max(x)-min(x),1);true(n_orderedP-max(x),1)]},orderedP_rc_u_loc(n_orderedP_rc_ui)')));
P_no_loop_mask=true(n_orderedP,1);
for loop=torow(orderedP_rc_u_loc(n_orderedP_rc_ui))
    lep = [min(loop{:}),max(loop{:})];
    if mod(lep(1)-lep(2),n_orderedP)<mod(lep(2)-lep(1),n_orderedP)
        P_no_loop_mask(1:(lep(1)-1))=false;
        P_no_loop_mask(lep(2):end)=false;
    else
        P_no_loop_mask((lep(1)+1):lep(2))=false;
    end
end
orderedP_rc_no_loop=orderedP_rc(P_no_loop_mask,:);

orderedP_i=sub2ind(globalsz,orderedP_rc_no_loop(:,1),orderedP_rc_no_loop(:,2));
Pmask=false(globalsz);
Pmask(orderedP_i)=true;
P_info.P_i=find(Pmask);
if numel(P_info.P_i)~=numel(orderedP_i)
    error('ordered P indices are not unique');
end
[~,P_info.PMi]=ismember(P_info.P_i,orderedP_i);
[~,P_info.PSi]=sort(P_info.PMi,'ascend');
P_info.P_rc=orderedP_rc_no_loop(P_info.PMi,:);
P_info.PMi=P_info.PMi./numel(P_info.PMi);
P_info.globalsz=globalsz;
Pfilled_mask = poly2mask(P_info.P_rc(P_info.PSi,2),P_info.P_rc(P_info.PSi,1),P_info.globalsz(1),P_info.globalsz(2));
% P_info.poly=s_mask2poly(Pfilled_mask,options.reduce_poly_length_to);

q_g_info=s_update_g_info_rot(fr_info{q}.g_info,TR(q,3));

pq_ex_rc = [sp_ex;sq_ex];
pq_ex_i = sub2ind(globalsz,pq_ex_rc(:,1),pq_ex_rc(:,2));
pq_ex_vals = [p_ex_only_vals;fr_info{q}.ex_only_vals];
pq_ex_vals_rgb = [p_ex_only_vals_rgb;fr_info{q}.ex_only_vals_rgb];
pq_ex_bins = [p_g_info.ex_bins;q_g_info.ex_bins];
pq_ex_mag = [p_g_info.ex_mag;q_g_info.ex_mag];
[P_info.P_ex_i,set2i]=unique(pq_ex_i);
P_info.P_ex_rc = pq_ex_rc(set2i,:);
P_info.P_ex_vals = pq_ex_vals(set2i,:);
P_info.P_ex_vals_rgb = pq_ex_vals_rgb(set2i,:);
P_info.g_info.ex_bins = pq_ex_bins(set2i);
P_info.g_info.ex_mag = pq_ex_mag(set2i);

[P_pmask,P_ploc]=ismember(P_info.P_rc,sp,'rows');
[P_qmask,P_qloc]=ismember(P_info.P_rc,sq,'rows');
[P_Pexmask,P_Pexloc]=ismember(P_info.P_i,P_info.P_ex_i);

P_info.P_vals(P_Pexmask,:)=P_info.P_ex_vals(P_Pexloc(P_Pexmask),:);
P_info.P_vals(P_qmask,:)=fr_info{q}.c_vals(P_qloc(P_qmask),:);
P_info.P_vals(P_pmask,:)=p_c_vals(P_ploc(P_pmask),:);
P_info.P_vals_rgb(P_Pexmask,:)=P_info.P_ex_vals_rgb(P_Pexloc(P_Pexmask),:);
P_info.P_vals_rgb(P_qmask,:)=fr_info{q}.c_vals_rgb(P_qloc(P_qmask),:);
P_info.P_vals_rgb(P_pmask,:)=p_c_vals_rgb(P_ploc(P_pmask),:);
P_info.g_info.c_bins(P_Pexmask)=P_info.g_info.ex_bins(P_Pexloc(P_Pexmask));
P_info.g_info.c_bins(P_qmask)=q_g_info.c_bins(P_qloc(P_qmask));
P_info.g_info.c_bins(P_pmask)=p_g_info.c_bins(P_ploc(P_pmask));
P_info.g_info.c_mag(P_Pexmask)=P_info.g_info.ex_mag(P_Pexloc(P_Pexmask));
P_info.g_info.c_mag(P_qmask)=q_g_info.c_mag(P_qloc(P_qmask));
P_info.g_info.c_mag(P_pmask)=p_g_info.c_mag(P_ploc(P_pmask));
P_info.g_info.c_mag=tocolumn(P_info.g_info.c_mag);
P_info.g_info.c_bins=tocolumn(P_info.g_info.c_bins);
P_info.g_info.binwidth=p_g_info.binwidth;
P_info.g_info.nbins=p_g_info.nbins;
% P_info.P_gdir(P_Pexmask,:)=P_info.P_ex_gdir(P_Pexloc(P_Pexmask),:);
% P_info.P_gdir(P_qmask,:)=fr_info{q}.gdir(P_qloc(P_qmask),:);
% P_info.P_gdir(P_pmask,:)=p_c_gdir(P_ploc(P_pmask),:);
% P_info.P_g_w8s(P_qmask,:)=fr_info{q}.g_w8s(P_qloc(P_qmask),:);
% P_info.P_g_w8s(P_pmask,:)=p_c_gw8s(P_ploc(P_pmask),:);

%DEBUG:
DebugPadding=20;
Debug_masks{1}=s_imtrim(Pmask,Pmask,DebugPadding);
Debug_masks{2}=false(globalsz);
Debug_masks{2}(sub2ind(globalsz,sp(:,1),sp(:,2)))=true;
Debug_masks{2}(sub2ind(globalsz,sq(:,1),sq(:,2)))=true;
Debug_masks{2}=s_imtrim(Debug_masks{2},Pmask,DebugPadding);
Debug_masks{3} = inf(globalsz);
Debug_masks{3}(Pmask)=P_info.PMi;
Debug_masks{3}=s_imtrim(Debug_masks{3},Pmask,DebugPadding);
Debug_masks{4}=Pfilled_mask;
%         figure(1);subplot(1,2,1);imshow(Pmask);subplot(1,2,2);imshow(pqmask);figure(2);imshow(convert_distance_color(coloredPmask));colormap(jet);colorbar;
tmp=1;
%-----
end
function [abs_tforms,calc_pairs_tforms_mask] = get_abs_tforms(rel_tforms,abs_tforms,pair_mask,upT,upR,C,calc_pairs_tforms_mask)
    n_pairs = size(C,1);
    for ci=1:n_pairs
        if (sum(pair_mask(ci,:))==1 && ~calc_pairs_tforms_mask(ci))
            n_scores = size(rel_tforms{ci},1);
            p=C(ci,pair_mask(ci,:));
            pair_rel_tforms = rel_tforms{ci};
            if (find(~pair_mask(ci,:))==1) %q==pair(1)
                pair_rel_tforms = s_invtform(pair_rel_tforms);
            end
            rep_p_tfroms = repmat([upT(p,:),upR(p)],n_scores,1);
            abs_tforms{ci} = s_compose_tforms(rep_p_tfroms,pair_rel_tforms);
            calc_pairs_tforms_mask(ci) = true;
        end
    end
end
function [is_valid_tform] = disable_tforms_cspace(abs_tforms,is_valid_tform,TR,pair_mask,pairs,cspace,Q)
for ci=1:size(pairs,1)
    if (sum(pair_mask(ci,:))==1)
        q=pairs(ci,~pair_mask(ci,:));
        if any(Q==q) % q may not be in Q when more then one image is extracted, then q may actually be an already placed fragment in a previous image
            tforms2check = abs_tforms{ci}(is_valid_tform{ci},:);
            is_valid_tforms2check=cspace_tforms_validation(tforms2check,pairs,q,pair_mask,cspace,TR,ci);
            is_valid_tform{ci}(is_valid_tform{ci}) = is_valid_tforms2check;        
        end
    end
end
end
function [is_valid]=cspace_tforms_validation(tforms2check,pairs,q,Ppairs_mask,cspace,TR,ci)
    if ~exist('ci','var')
        ci=0;
    end
    [rot_vals,~,pos]=unique(tforms2check(:,3));
    n_rot=numel(rot_vals);
    is_valid = true(size(tforms2check,1),1);
    [~,qpairs_i] = find(pairs'==q);
    qpairs_i = qpairs_i(sum(Ppairs_mask(qpairs_i,:),2)==1);
    if ~isrow(qpairs_i)
        qpairs_i = qpairs_i';
    end
    for qp_i=qpairs_i
        if ci~=qp_i
            p=pairs(qp_i,Ppairs_mask(qp_i,:));
            for r_i = 1:n_rot
                tform_i=find(pos==r_i);
                Rq=rot_vals(r_i);
                if any(tforms2check(tform_i,3)~= Rq)
                    error('All rotations in this set should be equal, Rq=%0.3f',Rq);
                end
                CSP=get_cspace(Ppairs_mask,qp_i,TR(p,:),Rq,cspace);
                csppoly.x = int64(CSP(:,1));
                csppoly.y = int64(CSP(:,2));
                for j=-10:0.1:-9
                    cspace_inset=clipper(csppoly,j,2);
                    [~,maxi]=max(arrayfun(@(z) numel(z.x),cspace_inset));
    %                         fprintf('ci=%d,qp_i=%d,r_i=%d\n',ci,qp_i,r_i);
                    cspace_inset_poly = [cspace_inset(maxi).x,cspace_inset(maxi).y];
                    if ~isempty(cspace_inset_poly)
                        break;
                    end
                end
                if ~isempty(cspace_inset_poly)
                    is_in_cspoly = c_inpoly(tforms2check(tform_i,2:-1:1),cspace_inset_poly);
                    is_valid(tform_i) = is_valid(tform_i) & ~is_in_cspoly;
                else
                    warning('could not perform inset for cspace');
                    
                end
            end
        end
    end
end
function [out_candT,out_fval] = refine_q_tform( P_info,fr_info,fr_input_info,P,q,initPT,qT,pairs,cspace,options )
rot = (-options.rot_refine_window:0.5:options.rot_refine_window)'+qT(3);
tr_range=(-options.tr_refine_window:options.tr_refine_window)';
[x,y]=meshgrid(1:numel(tr_range),1:numel(tr_range));
x=x(:);y=y(:);
tr = bsxfun(@plus,qT(1:2),[tr_range(x),tr_range(y)]);
[t_i,r_i]=meshgrid(1:size(tr,1),1:numel(rot));t_i=t_i(:);r_i=r_i(:);
tforms2check=[tr(t_i,:),rot(r_i)];
Ppairs_mask = ismember(pairs,P);
is_valid=cspace_tforms_validation(tforms2check,pairs,q,Ppairs_mask,cspace,initPT);
options.small_ov_percentage = 0;
dis_info = s_dissPq( P_info,fr_info,fr_input_info,P,q,initPT,tforms2check(is_valid,:),[],options );
vtforms2check = dis_info.qT;
disP = dis_info.disP;
disq = dis_info.disq;
if isempty(vtforms2check)
    vtforms2check = qT;
    disP = 0;
    disq = 0;
end
n_vtforms2check = size(vtforms2check,1);
Pqdists = zeros(n_vtforms2check,1);
if options.refine_alpha ~= 0
    if isempty(P_info)
        if numel(P)~=1
            error('No P_info but |P|==%d',numel(P)); 
        end
        Ppoly = fr_input_info{P}.poly;
    else
        Ppoly = P_info.poly;
    end
    qpoly_z = fr_input_info{q}.poly;
    qpoly_ex_z = fr_input_info{q}.poly_ex;
    parfor ti=1:n_vtforms2check
        qpoly = s_shift_poly(qpoly_z,vtforms2check(ti,:));
        qpoly_ex = s_shift_poly(qpoly_ex_z,vtforms2check(ti,:));
        Pqdists(ti) = polygons_dist(qpoly,qpoly_ex,Ppoly);
    end
    % qpolys = arrayfun(@(Ti) s_shift_poly(fr_input_info{q}.poly,vtforms2check(Ti,:)),(1:size(vtforms2check,1))','UniformOutput',false);
    % qpolys_ex = arrayfun(@(Ti) s_shift_poly(fr_input_info{q}.poly_ex,vtforms2check(Ti,:)),(1:size(vtforms2check,1))','UniformOutput',false);
    % Pqdists = cellfun(@(qpoly,qpoly_ex) polygons_dist(qpoly,qpoly_ex,Ppoly),qpolys,qpolys_ex);
end
scores = disP + disq + options.refine_alpha.*Pqdists;
[out_fval,bestscore_i]=min(scores);
out_candT = vtforms2check(bestscore_i,:);
fprintf('best score props: disP=%0.4f, disq=%0.4f, a*f(D(P,q))=%0.4f, score=%0.4f\n',disP(bestscore_i),disq(bestscore_i),options.refine_alpha*Pqdists(bestscore_i),scores(bestscore_i));
%DEBUG: qpoly = s_shift_poly(qpoly_z,vtforms2check(bestscore_i,:));figure;plot([Ppoly(1,:),Ppoly(1,1)],[Ppoly(2,:),Ppoly(2,1)],[qpoly(1,:),qpoly(1,1)],[qpoly(2,:),qpoly(2,1)]);
%DEBUG: qpoly_ex = s_shift_poly(qpoly_ex_z,vtforms2check(bestscore_i,:));polygons_dist(qpoly,qpoly_ex,Ppoly);
end
function d=polygons_dist(pol1,pol1_ex,pol2)
TOL = 1.0e-12;
in=c_inpoly(pol2',pol1_ex');
if any(in)
    dists = p_poly_dist(pol2(1,in),pol2(2,in),pol1(1,:),pol1(2,:));
    if any(dists<-TOL)
        d=abs(min(dists));
    else
        d=log(1+max(dists));
    end
else
    d=inf;
end
%DEBUG: figure;plot([pol1(1,:),pol1(1,1)],[pol1(2,:),pol1(2,1)],[pol1_ex(1,:),pol1_ex(1,1)],[pol1_ex(2,:),pol1_ex(2,1)],[pol2(1,:),pol2(1,1)],[pol2(2,:),pol2(2,1)],pol2(1,in),pol2(2,in),'+k');
end
function [conf,confP,confq,bst_tforms,Q2Pmap]=reassemb_confidence(P_info,abs_tforms,is_valid_tforms,dis_info,P,Q,input,T,R,fr_info,multiq_tforms,multiqtf_p,options)
tic;
n_Q=numel(Q);
[Pi,Qi]=meshgrid(1:numel(P),1:n_Q);
Pi=Pi(:);Qi=Qi(:);PP=P(Pi);QQ=Q(Qi);
if isrow(PP)
    PP=PP';
end
if isrow(QQ)
    QQ=QQ';
end
conf_dis_info.n_pics=dis_info.n_pics+1; 
P_idx=conf_dis_info.n_pics;
[~,opairi]=ismember(sort([PP,QQ],2),sort(dis_info.pairs,2),'rows');
o2pairi=setdiff(1:dis_info.n_pairs,opairi);
dis1 = cell(dis_info.n_pairs,1);
dis2 = cell(dis_info.n_pairs,1);
tforms = cell(dis_info.n_pairs,1);
% p_ov_ti = cell(dis_info.n_pairs,1);
% p1_ov = cell(dis_info.n_pairs,1);
p1_ov_lim = cell(dis_info.n_pairs,1);
% p2_ov = cell(dis_info.n_pairs,1);
p2_ov_lim = cell(dis_info.n_pairs,1);

dis1(opairi)=dis_info.dis1(opairi);
dis2(opairi)=dis_info.dis2(opairi);
tforms(opairi)=arrayfun(@(x) {abs_tforms{x}},opairi);
% p_ov_ti(opairi)=dis_info.p_ov_ti(opairi);
% p1_ov(opairi)=dis_info.p1_ov(opairi);
p1_ov_lim(opairi)=dis_info.p1_ov_lim(opairi);
% p2_ov(opairi)=dis_info.p2_ov(opairi);
p2_ov_lim(opairi)=dis_info.p2_ov_lim(opairi);


dis1(o2pairi)=dis_info.dis1(o2pairi);
dis2(o2pairi)=dis_info.dis2(o2pairi);
tforms(o2pairi)=arrayfun(@(x) {abs_tforms{x}},o2pairi);
% p_ov_ti(o2pairi)=dis_info.p_ov_ti(o2pairi);
% p1_ov(o2pairi)=dis_info.p1_ov(o2pairi);
p1_ov_lim(o2pairi)=dis_info.p1_ov_lim(o2pairi);
% p2_ov(o2pairi)=dis_info.p2_ov(o2pairi);
p2_ov_lim(o2pairi)=dis_info.p2_ov_lim(o2pairi);

PQ_pairs=dis_info.pairs;
Pq_pairs=[repmat(P_idx,n_Q,1),Q'];
dis1_add = cell(n_Q,1);
dis2_add = cell(n_Q,1);
tforms_add = cell(n_Q,1);
% p_ov_ti_add = cell(n_Q,1);
% p1_ov_add = cell(n_Q,1);
% p2_ov_add = cell(n_Q,1);
p1_ov_lim_add = cell(n_Q,1);
p2_ov_lim_add = cell(n_Q,1);
globalsz=cell2mat(arrayfun(@(qi) {s_calclimits(input,fr_info,[T,R],P,Q(qi))},(1:n_Q)'));
options.globalsz=max(globalsz,[],1);

for qi=1:n_Q
    q=Q(qi);
    n_multi_tforms_before_call = size(multiq_tforms{qi},1);
%     [ disP_mult,disq_mult,~,~,Pq_ov_ti_mult,p_ov_P_mult,p_ov_q_mult,multiq_tforms{qi},multiqtf_p{qi}] =...
%         s_measure_dis( fr_info,input,P,q,[T,R],multiq_tforms{qi},multiqtf_p{qi},options );
    [multi_dis_info] = s_dissPq( P_info,fr_info,input,P,q,[T,R],multiq_tforms{qi},multiqtf_p{qi},options );
    multiq_tforms{qi}=multi_dis_info.qT;
    multiqtf_p{qi}=multi_dis_info.qT_p;
    dis1_add{qi}=multi_dis_info.disP;
    dis2_add{qi}=multi_dis_info.disq;
    tforms_add{qi}=multiq_tforms{qi};
%     p_ov_ti_add{qi}=multi_dis_info.p_ov_ti;
%     p1_ov_add{qi}=multi_dis_info.p_ov_P;
%     p2_ov_add{qi}=multi_dis_info.p_ov_q;
    p1_ov_lim_add{qi}=multi_dis_info.p_ov_lim_P;
    p2_ov_lim_add{qi}=multi_dis_info.p_ov_lim_q;
%     dis1_add{qi}=disP_mult;
%     dis2_add{qi}=disq_mult;
%     tforms_add{qi}=multiq_tforms{qi};
%     p_ov_ti_add{qi}=Pq_ov_ti_mult;
%     p1_ov_add{qi}=p_ov_P_mult;
%     p2_ov_add{qi}=p_ov_q_mult;
    fprintf('after dis re-calc on frag cand %d for multi frag ov, saving %d/%d tforms for confidence\n',q,size(multiq_tforms{qi},1),n_multi_tforms_before_call);
end
conf_dis_info.dis1=[dis1;dis1_add];
conf_dis_info.dis2=[dis2;dis2_add];
conf_dis_info.tforms=[tforms;tforms_add];
% conf_dis_info.p_ov_ti=[p_ov_ti;p_ov_ti_add];
% conf_dis_info.p1_ov=[p1_ov;p1_ov_add];
% conf_dis_info.p2_ov=[p2_ov;p2_ov_add];
conf_dis_info.p1_ov_lim=[p1_ov_lim;p1_ov_lim_add];
conf_dis_info.p2_ov_lim=[p2_ov_lim;p2_ov_lim_add];
conf_dis_info.pairs=[PQ_pairs;Pq_pairs];
conf_dis_info.n_pairs=size(conf_dis_info.pairs,1);
conf_dis_info.valid_tform_mask=[is_valid_tforms';cellfun(@(x) {true(size(x))},dis1_add)];
fprintf('Preparing dis_info for confidence. time: %0.3f sec\n',toc);
tStartConf=tic;
% conf_options.closeness_threshold = options.closeness_threshold;
conf_options.iou_thresh = options.iou_thresh;
conf_options.gamma_L = options.gamma_L;
conf_options.gamma_H = options.gamma_H;
conf_options.save_top_conf_frag_num = options.save_top_conf_frag_num;
conf_options.selected_pairs = func_output(2,@ismember,sort([PP,QQ;Pq_pairs],2),sort(conf_dis_info.pairs,2),'rows');
% conf_options.selected_pairs = func_output(2,@ismember,sort([PP,QQ],2),sort(conf_dis_info.pairs,2),'rows');
conf_options.reltforms = dis_info.tforms;
input{P_idx}.rgb = options.curr_assemb;
input{P_idx}.mask = options.curr_assemb_mask;
[ ~,~,allconf,allmatch_info ] = confidence( conf_dis_info,conf_options,input );
fprintf('confidence calculation finished. time: %0.3f sec\n',toc(tStartConf));
max_n_results=50;
conf = -inf(n_Q,max_n_results);
confP = -inf(n_Q,max_n_results);
confq = -inf(n_Q,max_n_results);
Q2Pmap = cell(n_Q,max_n_results);
bst_tforms = zeros(n_Q,max_n_results,3);

for qi=1:n_Q
    allconftmp=allconf;
    qmask=conf_dis_info.pairs==Q(qi);
    q_loc_mask=any(qmask,2);
    allconftmp(~q_loc_mask,:)=-inf;
    allconftmp=allconftmp(:);
    [bestqconf,bestqconf_i]=sort(allconftmp,'descend');
    n_results = min(find(~isinf(bestqconf),1,'last'),max_n_results);
    conf(qi,1:n_results)=bestqconf(1:n_results);
    [bestqconf_i_pair,bestqconf_i_toplevel]=ind2sub(size(allconf),bestqconf_i);
    bestqpairs=conf_dis_info.pairs(bestqconf_i_pair(1:n_results),:)';
    confPq=cell2mat(arrayfun(@(pr,tl) {[allmatch_info{pr}.conf1(tl),allmatch_info{pr}.conf2(tl)]},tocolumn(bestqconf_i_pair(1:n_results)),tocolumn(bestqconf_i_toplevel(1:n_results))))';
    confP(qi,1:n_results)=confPq(bestqpairs~=Q(qi));
    confq(qi,1:n_results)=confPq(bestqpairs==Q(qi));
    Q2Pmap(qi,1:n_results)=num2cell(bestqpairs(bestqpairs~=Q(qi)));
    Q2Pmap(qi,1:n_results)=cellfun(@(op,pr,tl) {get_multip(op,multiqtf_p{qi},allmatch_info,pr,tl,P_idx)},Q2Pmap(qi,1:n_results),num2cell(bestqconf_i_pair(1:n_results)'),num2cell(bestqconf_i_toplevel(1:n_results)'));
    bst_tforms(qi,1:n_results,:)=cell2mat(arrayfun(@(pr,tl) {s_tform_mean([conf_dis_info.tforms{pr}(allmatch_info{pr}.tform_i1(tl),:);...
        conf_dis_info.tforms{pr}(allmatch_info{pr}.tform_i2(tl),:)])},bestqconf_i_pair(1:n_results),bestqconf_i_toplevel(1:n_results)));
end
end
function np=get_multip(op,multiqtf_p,allmatch_info,pr,tl,P_idx)
if op==P_idx
    np=unique([multiqtf_p{allmatch_info{pr}.tform_i1(tl)};multiqtf_p{allmatch_info{pr}.tform_i2(tl)}]);
else
    np=op;
end
end
function [top_i,qi,maxconf,qtform,bb_mask]=nextp(conf2,confP,confq,bst_tforms)
    epsilon=1e-14;
    bb_mask=confP>epsilon&confq>epsilon;%false(size(confP));%
    if false&&any(bb_mask(:))
        [maxconf,maxi]=max(conf2(bb_mask));
        [bb_qi,bb_topi]=find(bb_mask);
        qi=bb_qi(maxi);
        top_i=bb_topi(maxi);
    else
        top_i = 1;
        [maxconf,qi] = max(conf2(:,top_i));        
    end
    
    qtform = reshape(bst_tforms(qi,top_i,:),1,3);
end
function [tforms,tf_p]=calc_multi_frag_tform(pair_mask,pairs,P,Q,rot,sampling_res,TR, fullmask,masks,masks_ex,cspace,cspace_ex,is_valid_tforms)
to_clipper_struct=@(cs) struct('x',int64(cs(:,1)),'y',int64(cs(:,2)));
tic;
tforms = cell(numel(Q),1);
tf_p = cell(numel(Q),1);
max_masksz=max(cellfun(@(x) max(size(x)),masks));
finalsz = max(max_masksz,max(size(fullmask)));
sz = [finalsz,finalsz];
% fullmask = s_canvasSize(fullmask,sz);
Pbb_loc_rc = s_rotate_ind_loose(get_bb_loc(fullmask),size(fullmask),sz,0);
spairs = sort(pairs,2);
for qi=1:numel(Q)
    q=Q(qi);
%     qmask = s_canvasSize(masks{q},sz);
    qmask = masks{q};
    tf_p_qi = cell(numel(rot),1);
%     tmptf_p_qi = cell(numel(rot),1);
    tforms_qi = cell(numel(rot),1);
    
    [~,pi2qpairs_i] = ismember(sort([tocolumn(P),repmat(q,numel(P),1)],2),spairs,'rows');
    qpairs_cs = cspace(:,pi2qpairs_i);
    qpairs_cs_ex = cspace_ex(:,pi2qpairs_i);
    qpairs_mask = pair_mask(pi2qpairs_i,:);
    isvalidp=arrayfun(@(pairi) any(is_valid_tforms{pairi}),pi2qpairs_i);
    validpi=torow(find(isvalidp));
%     ticBytes(gcp);
    parfor ri=1:numel(rot)
%     for ri=1:numel(rot)
%         pqr_cs_fnl=[];
%         pqr_cs_ex_fnl=[];
        r=rot(ri);
%         qmask_rot=imrotate(qmask,r,'crop');
        qbb_loc_rc = s_rotate_ind_loose(get_bb_loc(qmask),size(qmask),sz,r);
%         possible_tr=calc_tr(fullmask,qmask_rot);
        possible_tr=calc_tr(Pbb_loc_rc,qbb_loc_rc);
        is_tr_ex_valid=false(size(possible_tr,1),numel(P));
        is_tr_physicaly_valid=true(size(possible_tr,1),numel(P));
        for pi=validpi
            p=P(pi);
%             [~,qp_i] = ismember(sort([p,q]),spairs,'rows');
%             pqr_cs=get_cspace(pair_mask, qp_i, TR(p,:), r,cspace);
            pqr_cs=get_cspace(qpairs_mask, pi, TR(p,:), r,qpairs_cs);
            for j=-5:0.1:-4
                pqr_cs_ins=clipper(to_clipper_struct(pqr_cs),j,2);
                [~,maxi]=max(arrayfun(@(z) numel(z.x),pqr_cs_ins));
                pqr_cs_fnl = [pqr_cs_ins(maxi).x,pqr_cs_ins(maxi).y];
                if ~isempty(pqr_cs_fnl)
                    break;
                end
            end
%             pqr_cs_ex=get_cspace(pair_mask, qp_i, TR(p,:), r,cspace_ex);
            pqr_cs_ex=get_cspace(qpairs_mask, pi, TR(p,:), r,qpairs_cs_ex);
            for j=-2:0.1:-1
                pqr_cs_ex_ins=clipper(to_clipper_struct(pqr_cs_ex),j,2);
                [~,maxi]=max(arrayfun(@(z) numel(z.x),pqr_cs_ex_ins));
                pqr_cs_ex_fnl = [pqr_cs_ex_ins(maxi).x,pqr_cs_ex_ins(maxi).y];
                if ~isempty(pqr_cs_ex_fnl)
                    break;
                end
            end
            n_cs_ex = size(pqr_cs_ex_fnl,1);
            n_cs = size(pqr_cs_fnl,1);

            c_cs_ex = [(1:n_cs_ex-1)', (2:n_cs_ex)'; n_cs_ex, 1];%n_cs_ex x 2
            c_cs = [(1:n_cs-1)', (2:n_cs)'; n_cs, 1];%n_cs x 2

            if isempty(pqr_cs_fnl) || isempty(pqr_cs_ex_fnl)
                error('empty cspace after inset polygon.');
            else
                is_tr_physicaly_valid(:,pi) = ~c_inpoly([possible_tr(:,2),possible_tr(:,1)],pqr_cs_fnl,c_cs);
                is_tr_ex_valid(:,pi) = c_inpoly([possible_tr(:,2),possible_tr(:,1)],pqr_cs_ex_fnl,c_cs_ex);
            end

%             figure(2);tmp=possible_tr(174164,:);plot(tmp(2),tmp(1),'.r',[mft_cs{p,q,ri}.cs_ex(:,1);mft_cs{p,q,ri}.cs_ex(1,1)],[mft_cs{p,q,ri}.cs_ex(:,2);mft_cs{p,q,ri}.cs_ex(1,2)],'b',[pqr_cs_fnl(:,1);pqr_cs_fnl(1,1)],[pqr_cs_fnl(:,2);pqr_cs_fnl(1,2)],'k')
%             figure(2);tmp=possible_tr(174164,:);plot(tmp(2),tmp(1),'.r',[pqr_cs_ex_fnl(:,1);pqr_cs_ex_fnl(1,1)],[pqr_cs_ex_fnl(:,2);pqr_cs_ex_fnl(1,2)],'b',[pqr_cs_fnl(:,1);pqr_cs_fnl(1,1)],[pqr_cs_fnl(:,2);pqr_cs_fnl(1,2)],'k')
%             figure(3);tmpi=174164;tmp=s_imcomb(double(fullmask),double(masks_ex{p})/2,fullmask,masks_ex{p},[T(p,:),R(p)],'assrc');tmp=tmp/2+s_imshift(qmask_rot,possible_tr(tmpi,:))/3;imshow(tmp);
        end
        multi_tr_mask = all(is_tr_physicaly_valid,2)&sum(is_tr_ex_valid,2)>1;
        [tr_valid_i,tr_valid_p]=find(is_tr_ex_valid(multi_tr_mask,:));
        
%         allvalidtr_p = cellfun(@(x) {P(x)'}, num2cell(is_tr_ex_valid(multi_tr_mask,:),2));
        
        if sum(multi_tr_mask)>20
            [valid_tr,N,M,rdata_i,cdata_i] = s_uniform_sampling(double(possible_tr(multi_tr_mask,:)),sampling_res,'new');

            data_i = sub2ind([N,M],rdata_i,cdata_i);
            
%             subs = cell2mat(cellfun(@(di,ps) {repmat(di,numel(ps),1)},num2cell(data_i),allvalidtr_p));
%             vals = cell2mat(allvalidtr_p);
%             tmpvalid_p = accumarray(subs,vals,[],@(x) {x});
%             tmptf_p_qi{ri} = tmpvalid_p(~cellfun(@isempty,tmpvalid_p));
            
            valid_p = accumarray(data_i(tr_valid_i),P(tr_valid_p),[],@(x) {x});
            tf_p_qi{ri} = valid_p(unique(data_i));
            valid_tforms = [valid_tr,repmat(r,size(valid_tr,1),1)];
            tforms_qi{ri} = valid_tforms;
        else
%             tmptf_p_qi{ri} = allvalidtr_p;
            
            tf_p_qi{ri} = cellfun(@(x) {P(x)'}, num2cell(is_tr_ex_valid(multi_tr_mask,:),2));
            valid_tr=double(possible_tr(multi_tr_mask,:));
            valid_tforms = [valid_tr,repmat(r,size(valid_tr,1),1)];
            tforms_qi{ri} = valid_tforms;
        end
    end
%     tocBytes(gcp);
%     tmp_tfp = cellfun(@unique,cat(1,tmptf_p_qi{:}),'UniformOutput',false);
    tf_p{qi} = cat(1,tf_p_qi{:});
    tforms{qi} = cell2mat(tforms_qi);
    tf_p{qi} = cellfun(@unique,tf_p{qi},'UniformOutput',false);
%     if ~isequal(tmp_tfp,tf_p{qi})
%         error('new tfp method fail');
%     end
end
% tmpqi=1;tmpi=2000;figure(3);imshow(s_imcomb(fullmask,masks{Q(tmpqi)},fullmask,masks{Q(tmpqi)},tforms{tmpqi}(tmpi,:)));
fprintf('multi fragment matching translations validation finished. time: %0.3f s\n',toc);
end
function [tforms,tf_p]=recalc_multifr_tforms(pair_mask,pairs,P,Q,TR,cspace_ex,is_valid_tforms,abs_tforms)
to_clipper_struct=@(cs) struct('x',int64(cs(:,1)),'y',int64(cs(:,2)));
tic;
tforms = cell(numel(Q),1);
tf_p = cell(numel(Q),1);
spairs = sort(pairs,2);
for qi=1:numel(Q)
    q=Q(qi);
    
    [~,pi2qpairs_i] = ismember(sort([tocolumn(P),repmat(q,numel(P),1)],2),spairs,'rows');
    qpairs_cs_ex = cspace_ex(:,pi2qpairs_i);
    qpairs_mask = pair_mask(pi2qpairs_i,:);
    isvalidp=arrayfun(@(pairi) any(is_valid_tforms{pairi}),pi2qpairs_i);
    validpi=torow(find(isvalidp));
    tforms2check = cell2mat(tocolumn(cellfun(@(t,v) t(v,:),abs_tforms(pi2qpairs_i),is_valid_tforms(pi2qpairs_i),'UniformOutput',false)));
    [urot,~,tf2urot_i]=unique(tforms2check(:,3));
    n_urot=numel(urot);
    tf_p_qi = cell(n_urot,1);
    tforms_qi = cell(n_urot,1);
%     ticBytes(gcp);
    parfor ri=1:n_urot
        r=urot(ri);
        possible_tr=tforms2check(tf2urot_i==ri,:);
        is_tr_ex_valid=false(size(possible_tr,1),numel(P));
        for pi=validpi
            p=P(pi);
            pqr_cs_ex=get_cspace(qpairs_mask, pi, TR(p,:), r,qpairs_cs_ex);
            for j=-2:0.1:-1
                pqr_cs_ex_ins=clipper(to_clipper_struct(pqr_cs_ex),j,2);
                [~,maxi]=max(arrayfun(@(z) numel(z.x),pqr_cs_ex_ins));
                pqr_cs_ex_fnl = [pqr_cs_ex_ins(maxi).x,pqr_cs_ex_ins(maxi).y];
                if ~isempty(pqr_cs_ex_fnl)
                    break;
                end
            end
            n_cs_ex = size(pqr_cs_ex_fnl,1);
            c_cs_ex = [(1:n_cs_ex-1)', (2:n_cs_ex)'; n_cs_ex, 1];%n_cs_ex x 2
            if isempty(pqr_cs_ex_fnl)
                error('empty cspace after inset polygon.');
            else
                is_tr_ex_valid(:,pi) = c_inpoly([possible_tr(:,2),possible_tr(:,1)],pqr_cs_ex_fnl,c_cs_ex);
            end
        end
        multi_tr_mask = sum(is_tr_ex_valid,2)>1;
        tf_p_qi{ri} = arrayfun(@(x) {tocolumn(P(is_tr_ex_valid(x,:)))},find(multi_tr_mask));
        tforms_qi{ri} = double(possible_tr(multi_tr_mask,:));
    end
%     tocBytes(gcp);
    tf_p{qi} = cat(1,tf_p_qi{:});
    tforms{qi} = cell2mat(tforms_qi);
end
fprintf('multi fragment matching translations validation finished. time: %0.3f s\n',toc);
end
function [tr]=calc_tr(p_ex_mask,q_nex_mask)
% [r,c] = find(q_nex_mask);
% min_nf_r = min(r);max_nf_r = max(r);
% min_nf_c = min(c);max_nf_c = max(c);

if size(p_ex_mask,2)==2
    if size(q_nex_mask,2)~=2
        error('bad input size. |q_nex_mask|=%s',mat2str(size(q_nex_mask)));
    end
    qbb_rc=q_nex_mask;
    pbb_rc=p_ex_mask;
else
    qbb_rc = get_bb_loc(q_nex_mask);
    pbb_rc = get_bb_loc(p_ex_mask);
end

min_nf_r = min(qbb_rc(:,1));max_nf_r = max(qbb_rc(:,1));
min_nf_c = min(qbb_rc(:,2));max_nf_c = max(qbb_rc(:,2));

n_nf_r_el = max_nf_r-min_nf_r+1;
n_nf_c_el = max_nf_c-min_nf_c+1;

% [r,c] = find(p_ex_mask);
% min_f_r = min(r);max_f_r = max(r);
% min_f_c = min(c);max_f_c = max(c);

min_f_r = min(pbb_rc(:,1));max_f_r = max(pbb_rc(:,1));
min_f_c = min(pbb_rc(:,2));max_f_c = max(pbb_rc(:,2));

lowest_start_r = (min_f_r - n_nf_r_el + 1);
lowest_start_c = (min_f_c - n_nf_c_el + 1);
dr=(lowest_start_r:max_f_r)';
dc=(lowest_start_c:max_f_c);
M = length(dr); N = length(dc);
dr=repmat(dr,1,N);
dc=repmat(dc,M,1);
stnfr = repmat(min_nf_r,M,N);
stnfc = repmat(min_nf_c,M,N);
set1 = dr-stnfr;
set2 = dc-stnfc;
tr(:,1) = set1(:);
tr(:,2) = set2(:);
end
function [bb_loc_rc] = get_bb_loc(mask)
[r,c] = find(mask);
min_r = min(r);max_r = max(r);
min_c = min(c);max_c = max(c);
bb_loc_rc=[min_r,min_c;
           min_r,max_c;
           max_r,min_c;
           max_r,max_c];
end
function [CSP]=get_cspace(pair_mask, qp_i, TRp, Rq,cspace)
Rp=TRp(3);
q_was_fixed_in_cs_calc = find(~pair_mask(qp_i,:))==1;%true -> q==pair(1) , false -> p == pair(1)
if q_was_fixed_in_cs_calc
    csangle = Rp-Rq;
    reverse_rotangle = 180-Rq;
else 
    csangle = -Rp+Rq;
    reverse_rotangle = -Rp;
end
csangle = mod(round(csangle),360);
reverse_rotangle = mod(reverse_rotangle,360);

reverseRot=rotz(reverse_rotangle);
reverseRot=reverseRot(1:2,1:2);

pqCSPoly = (reverseRot*cspace{csangle+1,qp_i})';
CSP = bsxfun(@plus,TRp(2:-1:1),pqCSPoly);
end

