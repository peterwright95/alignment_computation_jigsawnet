function [ globalsz ] = s_calclimits( fr_input_info,fr_add_info,initPT,P,q )
    limits = s_find_transl_limits(fr_input_info(P),initPT(P,:));
    diagP_l     = norm(limits);
    diagq_l    = norm(fr_input_info{q}.origsz_ex);
    sz          = max(round(diagP_l + 2*diagq_l),round(2*diagP_l + diagq_l));
    
%     max_masks_sz = max(s_calc_globalsz(fr_input_info,[tocolumn(P),repmat(q,numel(P))]));
%     max_masks_sz = max(arrayfun(@(p) max(size(fr_add_info{p}.c_mask)),[P,q]));
    finalsz = sz;%max(sz, max_masks_sz);

    globalsz = [finalsz,finalsz];
end

