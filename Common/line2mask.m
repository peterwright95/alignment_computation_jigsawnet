function [ xv,yv,mask ] = line2mask( x1,y1,x2,y2,sz )
% Distance (in pixels) between the two endpoints
nPoints = ceil(sqrt((x2 - x1).^2 + (y2 - y1).^2)) + 1;

% Determine x and y locations along the line
xv = round(linspace(x1, x2, nPoints))';
yv = round(linspace(y1, y2, nPoints))';
xy=unique([xv,yv],'rows','stable');
xv=xy(:,1);
yv=xy(:,2);
if exist('sz','var')
    % Replace the relevant values within the mask
    if numel(sz)==1
        sz=[sz,sz];
    end
    mask=false(sz);
    mask(sub2ind(sz, yv, xv)) = true;
else
    mask=[];
end
end

