function [cont_mask]=s_mask_contour(mask)
    mconv = [1,1,1;1,-8,1;1,1,1];
%     [cont_rc(:,1),cont_rc(:,2)] = find(conv2(double(~mask),mconv,'same')>0);
    cont_mask = conv2(double(~mask),mconv,'same')>0;
end

