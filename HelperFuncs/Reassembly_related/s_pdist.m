function [ D ] = s_pdist( p1,p2,mode,dim,w )
dist_t='cityblock';% 'euclidean', 'cityblock'
if ~exist('w','var')
    w=[1,1,1];
end
if numel(w)~=3 || any(w<0|w>1)
    error('Wrong input : w (weights)');
end
if iscolumn(w)
    w=w';
end
if strcmp(dist_t,'euclidean')
    w=sqrt(w);
end
N1=size(p1,1);
N2=size(p2,1);
W1=repmat(w,N1,1);
W2=repmat(w,N2,1);
if strcmp(mode, 'pairwise')
    D=pdist2(p1.*W1,p2.*W2,dist_t);
elseif strcmp(mode, '1to1')
    if ~isequal(size(p1),size(p2))
        error('size of p1 and p2 must be equal');
    end
    if ~exist('dim','var')
        dim=2;
    end
    if strcmp(dist_t,'euclidean')
        D = sqrt(sum((p1.*W1 - p2.*W2).^2,dim));
    elseif strcmp(dist_t,'cityblock')
        D = sum(abs(p1.*W1 - p2.*W2),dim);
    else
        error('unrecognized dist type');
    end
elseif strcmp(mode,'grad')
    gdistfunc = @(zi,zj) abs(sum(repmat(zi,size(zj,1),1).*zj,2));
%     gdistfunc = @(zi,zj) abs(sum(bsxfun(@times,zi,zj),2));
    D = pdist2(p1,p2,gdistfunc);
    % D = abs(1-pdist2(p1,p2,'cosine'));
    D=D/max(D(:));
else
    error('unrecognized mode');
end
end

