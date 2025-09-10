function [ coeff, err, is_y_parallel,xfit,yfit  ] = points2linecoeff( x,y)
is_y_parallel=false;
x=tocolumn(x);
y=tocolumn(y);
X=[ones(length(x),1),x];
if rank(X)==1
    [x,y]=deal(y,x);
    X=[ones(length(x),1),x];
    is_y_parallel=true;
end
b=X\y;
% if strcmp(mode,'linear')
%    b=X\y;
% elseif strcmp(mode,'svd')
%    [U,S,V]=svd(X);
%    s=diag(S);
%    b=(1/s(1)) * (U(:,1)'*y) * (V(:,1)')+...
%      (1/s(2)) * (U(:,2)'*y) * (V(:,2)');   
% else
%     error('unknown mode');
% end
% err = sum((y - yfitted).^2)/(sum((y - mean(y)).^2)+eps);
errs = abs(b(2).*x-y+b(1))./sqrt(b(2).^2+1);
err=mean(errs);
coeff = b;

xfit = x;
yfit = X*b;
if is_y_parallel
    [xfit,yfit]=deal(yfit,xfit);
end
end

