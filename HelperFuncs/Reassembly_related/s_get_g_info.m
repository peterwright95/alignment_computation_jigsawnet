function g_info = s_get_g_info(in_p,p_c_mask,p_ex_only_mask,nbins,r)
if ~exist('r','var')
    r=0;
end
if ~isequal(size(p_c_mask), size(p_ex_only_mask))
    error('masks should be the same size. |c_mask|=%s, |ex_mask|=%s', mat2str(size(p_c_mask)), mat2str(size(p_ex_only_mask)));
end
sz = size(p_c_mask);
binwidth=360/nbins;
grad = s_canvasSize(in_p.grad,sz);
gx = s_canvasSize(in_p.gx,sz);
gy = s_canvasSize(in_p.gy,sz);
% if r~=0
%     grad = imrotate(grad,r,'crop');
%     gx = imrotate(gx,r,'crop');
%     gy = imrotate(gy,r,'crop');
%     p_c = imrotate(p_c,r,'crop');
% end
g_info.c_mag = grad(p_c_mask);
g_c_ang = mod(atan2d(gy(p_c_mask),gx(p_c_mask))+r,360);
g_c_ang(g_c_ang==360)=0;
g_info.c_bins = floor(g_c_ang./binwidth)+1;

grad_ex = s_canvasSize(in_p.grad_ex,sz);
gx_ex = s_canvasSize(in_p.gx_ex,sz);
gy_ex = s_canvasSize(in_p.gy_ex,sz);
% if r~=0
%     grad_ex = imrotate(grad_ex,r,'crop');
%     gx_ex = imrotate(gx_ex,r,'crop');
%     gy_ex = imrotate(gy_ex,r,'crop');
%     p_ex_only = imrotate(p_ex_only,r,'crop');
% end
g_info.ex_mag = grad_ex(p_ex_only_mask);
g_ex_ang = mod(atan2d(gy_ex(p_ex_only_mask),gx_ex(p_ex_only_mask))+r,360);
g_ex_ang(g_ex_ang==360)=0;
g_info.ex_bins = floor(g_ex_ang./binwidth)+1;
g_info.nbins = nbins;
g_info.binwidth = binwidth;
end
