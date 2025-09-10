function [ P,noncentP ] = s_mask2poly( mask,reduce_poly_length_to )
    %S_MASK2POLY: return a closed polygon from mask reduced by a factor of
    %             reduce_poly_length_to [0-1] 

    if ~exist('reduce_poly_length_to','var')
        reduce_poly_length_to=0.2;
    end
    poly=mask2poly(mask);
%     if (numel(poly)~=1)
%         warning('fragment polygon has holes');
%     end
    [~,longest_i]=max(arrayfun(@(x) x.Length,poly));
    P = [poly(longest_i).X;poly(longest_i).Y];
%     tic;
%     [mask_c(:,2),mask_c(:,1)] = find(s_mask_contour(mask));
%     [c_order_d,c_order]=pdist2(P',mask_c,'euclidean','Smallest',1);
%     c_order_d = c_order_d./max(c_order_d);
%     c_order_d = c_order_d + c_order;
%     [~,final_order] = sort(c_order_d);
%     toc;
    P = reduce_poly(P, round(poly(longest_i).Length*reduce_poly_length_to));
    [P(1,:),P(2,:)] = poly2ccw(P(1,:),P(2,:));
    noncentP=P;
    center = flip(size(mask)/2)'; % flipped because rows = y axis, cols = x axis
%     P=P-repmat(center,1,size(P,2));
    P=bsxfun(@minus,P,center);
    P=[P(:,1),unique(P(:,2:end)','rows','stable')']; % eliminate closed loops or repeating vertices
    
    %Debug: Test for repeating vertices: find(sum(abs(diff(tmp,1,2)),1)==0)
%     P(:,end)=[];
end

