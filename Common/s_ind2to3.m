function [ ind3 ] = s_ind2to3( ind2_, sz )
ind2=ind2_;
if (nargin==1)
    if (~islogical(ind2))
        error('if sz is not given then ind2 should be a mask');
    end
    sz = size(ind2);
    ind2 = find(ind2_);
end
if (size(sz)<2)
    error('wrong number of dimensions for ind2to3');
end
if isrow(ind2)
    ind2 = ind2';
end
szprod = sz(1)*sz(2);
% szprod_rep = repmat(szprod,numel(ind2),1);
% ind3=cumsum([ind2, szprod_rep ,szprod_rep],2);
ind3=[ind2,bsxfun(@plus,ind2,szprod),bsxfun(@plus,ind2,2*szprod)];
end

