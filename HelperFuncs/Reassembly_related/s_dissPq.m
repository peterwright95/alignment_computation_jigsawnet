function [dis_info,Pdata,qdata] = s_dissPq( varargin )
%   call 1 (recalculate data): [dis_info]=s_dissPq(P_info,fr_info,fr_input_info,P,q,initPT,qT,qT_p,options);
%   call 2 (use previously calculated data): [dis_info]=s_dissPq(Pdata,qdata,qT,qT_p,options);
%
%

    [Pdata,qdata,qT,qT_p,options]=prepare_data(varargin);
    if ~isequal(Pdata.globalsz,qdata.globalsz)
        error('different globalsz detected. Pdata.globalsz=%s,qdata.globalsz=%s',mat2str(Pdata.globalsz),mat2str(qdata.globalsz));
    else
        globalsz = Pdata.globalsz;
    end
%     indmask=reshape(1:prod(globalsz),globalsz);
    n_tforms = size(qT,1);
    disP = zeros(n_tforms,1);
    disq = zeros(n_tforms,1);
    ov_limitsP = nan(n_tforms,2);
    ov_limitsq = nan(n_tforms,2);
    n_ovP = zeros(n_tforms,1);
    n_ovq = zeros(n_tforms,1);

    [urot,~,rloc]=unique(mod(qT(:,3),360));
    n_urot=numel(urot);

    
    small_ov_thresh = min(numel(Pdata.Mi),numel(qdata.Mi)) * options.small_ov_percentage;
    for ri=1:n_urot
        r=urot(ri);
        ri_T_loc = find(rloc==ri);
        curr_rot_T = qT(ri_T_loc,:);
        inv_curr_rot_T = s_invtform(curr_rot_T);

        Prot_c_rc = s_rotate_ind_loose( Pdata.c_rc,Pdata.orig_sz,globalsz,-r );
        qrot_c_rc = s_rotate_ind_loose( qdata.c_rc,qdata.orig_sz,globalsz,r );
        
%         Prot_c_rc = s_rotate_ind( Pdata.c_i,indmask,globalsz,-r );
%         qrot_c_rc = s_rotate_ind( qdata.c_i,indmask,globalsz,r );
        
        P_r_g_info = s_update_g_info_rot(Pdata.g_info,-r);
        q_r_g_info = s_update_g_info_rot(qdata.g_info,r);
        
        [ curr_disP,~,n_ovP(ri_T_loc),~,curr_ov_limP ] = s_dissimilarity( qdata.ex_only_vals,Pdata.c_vals,qdata.ex_i,Prot_c_rc,[],inv_curr_rot_T(:,1:2),Pdata.no_ov_val_P,globalsz,Pdata.Mi,qdata.g_info,P_r_g_info,small_ov_thresh,options.smoothing_deg_window);
        [ curr_disq,~,n_ovq(ri_T_loc),~,curr_ov_limq ] = s_dissimilarity( Pdata.ex_only_vals,qdata.c_vals,Pdata.ex_i,qrot_c_rc,[],curr_rot_T(:,1:2),qdata.no_ov_val_q,globalsz,qdata.Mi,Pdata.g_info,q_r_g_info,small_ov_thresh,options.smoothing_deg_window);

        ov_limitsP(ri_T_loc,:) = curr_ov_limP;
        ov_limitsq(ri_T_loc,:) = curr_ov_limq;
        disP(ri_T_loc)=curr_disP(:,options.dis_version);
        disq(ri_T_loc)=curr_disq(:,options.dis_version);
    end
    legal_tforms_mask = ~isinf(disP)&~isinf(disq);
    dis_info.p_ov_lim_P=ov_limitsP(legal_tforms_mask,:);
    dis_info.p_ov_lim_q=ov_limitsq(legal_tforms_mask,:);
    dis_info.disP = disP(legal_tforms_mask);
    dis_info.disq = disq(legal_tforms_mask);
    dis_info.n_ovP = n_ovP(legal_tforms_mask);
    dis_info.n_ovq = n_ovq(legal_tforms_mask);
    dis_info.qT=qT(legal_tforms_mask,:);
    if ~isempty(qT_p)
        dis_info.qT_p=qT_p(legal_tforms_mask,:);
    end
end
function [Pdata,qdata,qT,qT_p,options]=prepare_data(arguments)
    if numel(arguments)==5
        % [dis_info]=s_dissPq(Pdata,qdata,qT,qT_p,options);
        Pdata=arguments{1};
        qdata=arguments{2};
        qT=arguments{3};
        qT_p=arguments{4};
        options=arguments{5};
    else
        % [dis_info]=s_dissPq(P_info,fr_info,fr_input_info,P,q,initPT,qT,qT_p,options);
        P_info=arguments{1};
        fr_info=arguments{2};
        fr_input_info=arguments{3};
        P=arguments{4};
        q=arguments{5};
        initPT=arguments{6};
        qT=arguments{7};
        qT_p=arguments{8};
        options=arguments{9};
        globalsz=s_calclimits(fr_input_info,fr_info,initPT,P,q);
        Pdata=prepare_Pdata(P_info,P,fr_input_info,fr_info,globalsz);
        qdata=prepare_frdata(fr_input_info,fr_info,q,Pdata.globalsz);
        non_ov_val = s_calc_non_ov_val(Pdata,qdata);
        Pdata.no_ov_val_P=non_ov_val(1);
        qdata.no_ov_val_q=non_ov_val(2);
    end    
end
function [Pdata]=prepare_Pdata(P_info,P,fr_input_info,fr_info,globalsz)
    Pdata=[];
    if ~isempty(P_info)
        c_rc=P_info.P_rc; %adjust to globalsz
        ex_rc=P_info.P_ex_rc; %adjust to globalsz
        if any(globalsz>P_info.globalsz)
           offset=round((globalsz-P_info.globalsz)./2);
           c_rc=bsxfun(@plus,c_rc,offset);
           ex_rc=bsxfun(@plus,ex_rc,offset);
           Pdata.globalsz=globalsz;
        else
           Pdata.globalsz=P_info.globalsz;
        end
        Pdata.c_rc=c_rc;
        Pdata.orig_sz = Pdata.globalsz;
        Pdata.c_i=sub2ind(Pdata.globalsz,c_rc(:,1),c_rc(:,2));
        Pdata.ex_i=sub2ind(Pdata.globalsz,ex_rc(:,1),ex_rc(:,2));
        Pdata.Mi=P_info.PMi;
        Pdata.ex_only_vals=P_info.P_ex_vals;
        Pdata.ex_only_vals_rgb=P_info.P_ex_vals_rgb;
%         Pdata.ex_gdir=P_info.P_ex_gdir;
        Pdata.c_vals=P_info.P_vals;
        Pdata.c_vals_rgb=P_info.P_vals_rgb;
        Pdata.g_info=P_info.g_info;
    else
        if numel(P)~=1
            error('P_info is empty but P contains more than 1 fragment');
        end
        [Pdata]=prepare_frdata(fr_input_info,fr_info,P,globalsz);
    end
end
function [data]=prepare_frdata(fr_input_info,fr_info,q,globalsz)
    data=[];
    data.c_rc = mask_indices(fr_info{q}.c_mask);
    data.orig_sz = size(fr_info{q}.c_mask);
%     data.c_i = find(s_canvasSize(fr_info{q}.c_mask,globalsz));
    [~,data.c_i] = mask_indices(fr_info{q}.c_mask,globalsz);
    [~,data.ex_i] = mask_indices(fr_info{q}.ex_only_mask,globalsz);
%     data.ex_i = find(s_canvasSize(fr_info{q}.ex_only_mask,globalsz));
    data.Mi=fr_input_info{q}.mask_Mi;
    data.ex_only_vals = fr_info{q}.ex_only_vals;
    data.ex_only_vals_rgb = fr_info{q}.ex_only_vals_rgb;
%     data.ex_gdir = fr_info{q}.gdir_ex;
    data.c_vals=fr_info{q}.c_vals;    
    data.c_vals_rgb=fr_info{q}.c_vals_rgb;
    data.g_info=fr_info{q}.g_info;
    data.globalsz = globalsz;
end