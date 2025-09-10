function [ ov_limits ] = s_ov2limits( ov,mask_Mi )
    if numel(mask_Mi)~=numel(ov)
        error('numel(mask order) (%d) ~= numel(overlap) (%d)',numel(mask_Mi),numel(ov));
    end
    ov_Mi=mask_Mi(ov);
    if isempty(ov_Mi)
        ov_limits=nan(1,2);
    else
        if any(ov_Mi>1|ov_Mi<0)
            error('overlap parameteric values must remain in [0,1] interval. instead was: [%s]',num2str(ov_Mi,' '));
        end
        ov_Mi_s=torow(sort(ov_Mi,'ascend'));
        ov_Mi_s_extended =[ov_Mi_s,ov_Mi_s(1)]; 
        seg_diff=cyclic_diff(ov_Mi_s_extended);
        [~,max_seg_diff_i]=max(seg_diff);
        ov_limits=flip(ov_Mi_s_extended(max_seg_diff_i:max_seg_diff_i+1));
    end
end

