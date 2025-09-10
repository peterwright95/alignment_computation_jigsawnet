function [ trim, bbox, trimmed_mask ] = s_imtrim( m,mask,padding )
if ~exist('padding','var')
    padding = 0;
end
[r,c] = find(mask);
bbox.xmin = min(c)-padding; 
bbox.ymin = min(r)-padding;
xmax = max(c); 
ymax = max(r);
bbox.width = xmax - bbox.xmin + 1 + padding; 
bbox.height = ymax - bbox.ymin + 1 + padding;
trim = imcrop(m, [bbox.xmin bbox.ymin bbox.width bbox.height]);
if (nargout>=3)
    trimmed_mask = imcrop(mask, [bbox.xmin bbox.ymin bbox.width bbox.height]);
end
end

