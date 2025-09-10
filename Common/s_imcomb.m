function [ com, com_mask,shifted_sec,shifted_sec_mask ] = s_imcomb( src, sec, src_mask, sec_mask, T,szmode)
    %S_IMCOMB Combine 2 images according to their contour
    if ~exist('szmode','var')
        szmode = 'fit';
    end
    isshifted_sec_req = nargout>2;
    shifted_sec = [];
    shifted_sec_mask = [];
    if (~exist('T','var'))
        T = [0,0,0];
    end
    if ~isrow(T)
        T = T';
    end
    if (any(T~=0))
        sec = imrotate(sec,T(3));
        sec_mask = imrotate(sec_mask,T(3));
        if strcmp(szmode,'fit')
            finalsz = max(size(src_mask),size(sec_mask));
            cent = size(sec_mask)/2;
            added_sz = T(1:2) + size(sec_mask) - cent + finalsz/2;
            beyond_wh_i = find(added_sz > finalsz);
            if ~isempty(beyond_wh_i)
                finalsz(beyond_wh_i) = finalsz(beyond_wh_i) + 2*(added_sz(beyond_wh_i)-finalsz(beyond_wh_i));
            end
            added_sz = T(1:2) + [1,1] - cent + finalsz/2;
            beyond_0_i = find(added_sz(1:2) < 1);
            if ~isempty(beyond_0_i)
                finalsz(beyond_0_i) = finalsz(beyond_0_i) + 2*(1 - added_sz(beyond_0_i));
            end
        end
    else
        finalsz = max(size(src_mask),size(sec_mask));
    end
    if strcmp(szmode,'fit')
        finalsz = round(finalsz);
    elseif strcmp(szmode,'assrc')
        finalsz = size(src_mask);
    else
        error('unknown size mode');
    end
    if isempty(src)
        com = sec;
        com_mask = sec_mask;
    else
        src=s_canvasSize(src,finalsz);
        src_mask=s_canvasSize(src_mask,finalsz);
        com = src;
        com_mask = src_mask;
        ind2d = find(sec_mask);
%         [ind2d_rc(:,1),ind2d_rc(:,2)] = find(sec_mask);
%         sind2d_rc = round(ind2d_rc ...
%                     + repmat(T(1:2),size(ind2d_rc,1),1)...
%                     - repmat(ceil(size(sec_mask)/2),size(ind2d_rc,1),1)...
%                     + repmat(ceil(finalsz/2),size(ind2d_rc,1),1));
        sind2d_rc = round(bsxfun(@plus,mask_indices(sec_mask,finalsz),T(1:2)));
%         sind2d_rc = ceil(bsxfun(@plus,ind2d_rc,T(1:2)-size(sec_mask)/2+finalsz/2));
        if any(sind2d_rc(:,1)>finalsz(1))||...
                any(sind2d_rc(:,2)>finalsz(2))||...
                any(sind2d_rc(:,1)<1)||...
                any(sind2d_rc(:,2)<1)
            error('shifted indices out of bounds');
        end
        sind2d_i = sub2ind(finalsz, sind2d_rc(:,1),sind2d_rc(:,2));
        if size(com,3)==1
            ind = ind2d;
            com_ind = sind2d_i;
        else
            ind = s_ind2to3(sec_mask);
            com_ind = s_ind2to3(sind2d_i, finalsz);
        end
        com(com_ind)=sec(ind);
        com_mask(sind2d_i)=sec_mask(ind2d);
        if isshifted_sec_req 
            shifted_sec = zeros(size(com),'uint8');
            shifted_sec(com_ind) = sec(ind);
            shifted_sec_mask = false(size(com_mask));
            shifted_sec_mask(sind2d_i) = true;
        end
        
        if size(com,3)~=1
            ind = s_ind2to3(~com_mask);
            com(ind) = 255;
        end
        
    end
end

