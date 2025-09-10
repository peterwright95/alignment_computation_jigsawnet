function [dis_info]=s_calcerrors( fr_input,fr_info,options)
% Init workspace
rot = options.rot;
sampling_res = options.sampling_res;
smoothing_deg_window = options.smoothing_deg_window;
n_rot = numel(rot);
N=numel(fr_input);
pairs = nchoosek(1:N, 2);
n_pairs = size(pairs,1);
dis_info.dis_ov1 = cell(1, n_pairs);
dis_info.dis_ov2 = cell(1, n_pairs);
dis_info.dis1 = cell(1, n_pairs);
dis_info.dis2 = cell(1, n_pairs);
dis_info.tforms = cell(1, n_pairs);
dis_info.n_ov1 = cell(1, n_pairs);
dis_info.n_ov2 = cell(1, n_pairs);
dis_info.non_ov_val = zeros(N,N);
dis_info.cspace = cell(n_rot,n_pairs);
dis_info.cspace_ex = cell(n_rot,n_pairs);
dis_info.p_ov_ti=cell(1, n_pairs);
dis_info.p1_ov=cell(1, n_pairs);
dis_info.p2_ov=cell(1, n_pairs);
dis_info.rot=rot;
dis_info.pairs = pairs;
globalsz = s_calc_globalsz(fr_input,pairs);
[allcspace,allcspace_ex] = s_calc_allcspace(fr_input,pairs);

for ci=1:n_pairs

    p1=pairs(ci,1);
    p2=pairs(ci,2);
    in_p1=fr_input{p1};
    in_p2=fr_input{p2};
    tic;
    origsz1 = size(in_p1.mask);
    origsz2 = size(in_p2.mask);
    sz = [globalsz(ci),globalsz(ci)];
   
    b_p1_ext_z_rc = mask_indices(in_p1.mask_ex,sz);
    
    tforms      = cell(n_rot,1);
    dis1       = cell(n_rot,1);
    dis2       = cell(n_rot,1);
    n_ov1       = cell(n_rot,1);
    n_ov2       = cell(n_rot,1);
    ov_limits1 = cell(n_rot,1);
    ov_limits2 = cell(n_rot,1);
    orig_n_tforms=zeros(n_rot,1);
    n_valid_tforms=zeros(n_rot,1);
        
    p1val=fr_info{p1}.c_vals;
    p1val_ex=fr_info{p1}.ex_only_vals;
    p1_g_info=fr_info{p1}.g_info;
    p2val=fr_info{p2}.c_vals;
    p2val_ex=fr_info{p2}.ex_only_vals;
    p2_g_info=fr_info{p2}.g_info;
    
    no_ov_val = s_calc_non_ov_val(fr_info{p1}, fr_info{p2});
    no_ov_val_on_1=no_ov_val(1);
    no_ov_val_on_2=no_ov_val(2);
    p1_z_c_rc = mask_indices(fr_info{p1}.c_mask);
    p2_z_c_rc = mask_indices(fr_info{p2}.c_mask);

    [~,p1_z_ex_only_i] = mask_indices(fr_info{p1}.ex_only_mask,sz);
    [~,p2_z_ex_only_i] = mask_indices(fr_info{p2}.ex_only_mask,sz);

    p1curve = in_p1.mask_Mi;
    p2curve = in_p2.mask_Mi;
    ci_cspace = allcspace(round(rot)+1,ci);
    ci_cspace_ex = allcspace_ex(round(rot)+1,ci);
    small_ov_thresh = min(numel(p1curve),numel(p2curve)) * options.small_ov_percentage;

% ticBytes(gcp);
    parfor ri = 1:n_rot 
%     for ri = 1:n_rot 
        r=rot(ri);

        p1_c_rc = s_rotate_ind_loose( p1_z_c_rc,origsz1,sz,-r );
        p2_c_rc = s_rotate_ind_loose( p2_z_c_rc,origsz2,sz,r );

        K = calc_translations(b_p1_ext_z_rc, p2_c_rc,ci_cspace{ri},ci_cspace_ex{ri});

        KIcurRot=s_uniform_sampling(double(K),sampling_res,'new');
                
        curr_rot_tforms = double([KIcurRot,repmat(r,size(KIcurRot,1),1)]); 
        curr_invtforms = s_invtform(curr_rot_tforms);
        
        orig_n_tforms(ri) = size(curr_rot_tforms,1);
        
        p1_r_g_info=s_update_g_info_rot(p1_g_info,-r);
        p2_r_g_info=s_update_g_info_rot(p2_g_info,r);

        
        [ curr_dis1,~,curr_n_ov1,~,curr_ov_limits1] = s_dissimilarity( p2val_ex,p1val,p2_z_ex_only_i,p1_c_rc,[],curr_invtforms(:,1:2),no_ov_val_on_1,sz,p1curve,p2_g_info,p1_r_g_info,small_ov_thresh,smoothing_deg_window);
        [ curr_dis2,~,curr_n_ov2,~,curr_ov_limits2] = s_dissimilarity( p1val_ex,p2val,p1_z_ex_only_i,p2_c_rc,[],curr_rot_tforms(:,1:2),no_ov_val_on_2,sz,p2curve,p1_g_info,p2_r_g_info,small_ov_thresh,smoothing_deg_window);
        
        legal_tforms_mask = ~isinf(curr_dis1(:,1))&~isinf(curr_dis2(:,1));
        n_valid_tforms(ri)=nnz(legal_tforms_mask);
        dis1{ri} = curr_dis1(legal_tforms_mask,:);
        dis2{ri} = curr_dis2(legal_tforms_mask,:);
        ov_limits1{ri}=curr_ov_limits1(legal_tforms_mask,:);
        ov_limits2{ri}=curr_ov_limits2(legal_tforms_mask,:);
        tforms{ri}=curr_rot_tforms(legal_tforms_mask,:);
        n_ov1{ri}=curr_n_ov1(legal_tforms_mask,:);
        n_ov2{ri}=curr_n_ov2(legal_tforms_mask,:);
    end
% tocBytes(gcp);
    dis_info.tforms{ci}     = cell2mat(tforms);
    dis_info.dis1{ci}       = cell2mat(dis1);
    dis_info.dis2{ci}       = cell2mat(dis2);
    dis_info.n_ov1{ci}      = cell2mat(n_ov1);
    dis_info.n_ov2{ci}      = cell2mat(n_ov2);
    dis_info.cspace         = allcspace;
    dis_info.cspace_ex      = allcspace_ex;
    dis_info.non_ov_val(p1,p2)=no_ov_val_on_1;
    dis_info.non_ov_val(p2,p1)=no_ov_val_on_2;
    dis_info.p1_ov_lim{ci} = cell2mat(ov_limits1);
    dis_info.p2_ov_lim{ci} = cell2mat(ov_limits2);
    

    
    
    
    fprintf('Disimilarity calc progress: (%d/%d), %0.3f%%, calc-time: %0.3f sec\n',ci,n_pairs,ci/n_pairs*100, toc);
    fprintf('Saving %d tforms out of %d total tforms (valid: %d)\n', size(dis_info.tforms{ci},1),sum(orig_n_tforms), sum(n_valid_tforms));
end

end

function [tr]=calc_translations(f_ex_rc,nf_nex_rc,cspace,cspace_ex)
% [r,c] = find(nf_nex_mask);
r=nf_nex_rc(:,1);
c=nf_nex_rc(:,2);
min_nf_r = min(r);max_nf_r = max(r);
min_nf_c = min(c);max_nf_c = max(c);
n_nf_r_el = max_nf_r-min_nf_r+1;
n_nf_c_el = max_nf_c-min_nf_c+1;

% [r,c] = find(f_ex_mask);
r=f_ex_rc(:,1);
c=f_ex_rc(:,2);
min_f_r = min(r);max_f_r = max(r);
min_f_c = min(c);max_f_c = max(c);

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
set(:,1) = set1(:);
set(:,2) = set2(:);

n1 = size(cspace_ex,2);
n2 = size(cspace,2);

c1 = [(1:n1-1)', (2:n1)'; n1, 1];
c2 = [(1:n2-1)', (2:n2)'; n2, 1];

node  = [cspace_ex'; cspace'];
cnect = [c1;c2+n1];
    
is_in_cspace = c_inpoly([set(:,2),set(:,1)],node,cnect);
tr = set(is_in_cspace,:);
end
%-----------------------------