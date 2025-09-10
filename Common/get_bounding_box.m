function [ bb ] = get_bounding_box( mask )
[r,c]=find(mask);
bb.rmin=min(r);
bb.rmax=max(r);
bb.cmin=min(c);
bb.cmax=max(c);
end

