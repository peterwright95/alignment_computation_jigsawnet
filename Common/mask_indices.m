function [mask_rc,mask_i]=mask_indices(mask, outmasksz)
    [mask_rc(:,1),mask_rc(:,2)] = find(mask);
    insz = size(mask);
    if exist('outmasksz','var')
        incn = insz./2;
        outcn = outmasksz./2;
        mask_rc = bsxfun(@plus,mask_rc,ceil(outcn-incn));
    else
        outmasksz = insz;
    end
    if nargout>1
        mask_i = sub2ind(outmasksz,mask_rc(:,1),mask_rc(:,2));
    end
end