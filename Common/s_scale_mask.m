function [ nm ] = s_scale_mask( m, n_px, mode)
% SCALE_MASK - scale the input mask by a specific number of pixels
% input:
% m - the mask of the fragment (logical values: 1 | 0). has to include
%     exactly one connected component
% n_px - the number of pixels to move the contour in the direction of its
%        sign. (e.g. n_px=-2 means to return a mask smaller by 2 pixel in
%        the normal direction of each pixel.
% mode - 'dilate' , 'conv2'(default)
%
    
    % dilate shape options : 'arbitrary', 'square', 'diamond', 'rectangle',
    %                        'octagon', 'line', 'disk', 'sphere', 'cube', 'cuboid'
    dltshape = 'square';
    if ~exist('mode','var')
        mode = 'conv2';
    end

    if ~islogical(m)%(sum(sum(m~=0 & m~=1))>0)
        error('Invalid input image'); 
    end
    
    if n_px ~= 0
        abs_px = abs(n_px);
        sign = n_px/abs_px;

        if strcmp(mode,'dilate')
            if (sign>0)
                nm = imdilate(m, strel(dltshape, abs_px*2));
            else
                nm = ~imdilate(~m, strel(dltshape, abs_px*2));
            end
        elseif strcmp(mode,'conv2')
            out_val = m(1);
            in_val = 1-out_val;

            dmask = single(m);
            mconv = [1,1,1;1,-8,1;1,1,1];
            for i=1:abs_px
                if (sign>0)
                    dmask(conv2(dmask,mconv,'same')>0) = in_val;
                else
                    dmask(conv2(1-dmask,mconv,'same')>0) = out_val;
                end
            end
            nm = cast(dmask, class(m));
        else
            error('unknown mode for scale mask');
        end
    else
        nm=m;
    end
    
    [bmlbls,bmlbls_n]=bwlabel(nm,4);
    if bmlbls_n>1
        bmcomp_sizes=accumarray(bmlbls(bmlbls~=0),1);
        [~,maxlbls_i] = max(bmcomp_sizes);
%         fprintf('[s_scale_mask]-%d connected components of sizes %s detected, using component %d\n',bmlbls_n,mat2str(bmcomp_sizes),maxlbls_i);
        nm = bmlbls==maxlbls_i;
    end
end

