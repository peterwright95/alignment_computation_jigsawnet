function [ out_rc ] = s_rotate_ind_loose( in_rc,in_sz,out_sz,r )
mode = 'rotmat';
if strcmp(mode,'rotmat')
    in_cn = in_sz./2;
    out_cn = out_sz./2;
    out_rc = round((in_rc-in_cn)*rot2D(r)'+out_cn);
else
    after_rot_sz=torow(ceil(abs(rot2D(r))*in_sz'));
    loose_sz = max(after_rot_sz,in_sz);
    indmask = reshape(1:prod(loose_sz),loose_sz);
    rotindmask=imrotate(indmask,-r,'nearest','crop');
    in_cn = in_sz./2;
    loose_cn = loose_sz./2;
    loose_rc = bsxfun(@plus,in_rc,ceil(loose_cn-in_cn));
    loose_i = sub2ind(loose_sz,loose_rc(:,1),loose_rc(:,2));
    [rot_rc(:,1),rot_rc(:,2)]=ind2sub(loose_sz,rotindmask(loose_i));
    out_cn = out_sz./2;
    out_rc = bsxfun(@plus,rot_rc,ceil(out_cn-loose_cn));
end
end

