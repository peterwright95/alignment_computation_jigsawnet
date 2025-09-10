function [r,maxi]=infmax(d,ind)
notinfi=find(~isinf(d(ind)));
[r,maxi]=nanmax(d(ind(notinfi)));
maxi=ind(notinfi(maxi));
if isempty(r)
    r=-inf;
    maxi=0;
end
end