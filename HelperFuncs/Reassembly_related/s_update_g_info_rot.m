function [ g_info_r ] = s_update_g_info_rot(g_info,r)
bin_step = round(r/g_info.binwidth);
% if bin_step ~= floor(bin_step)
%     error('rotation should divide by binwidth, r=%d,width=%d, step=%0.4f',r,g_info.binwidth,bin_step);
% end
g_info_r=g_info;
g_info_r.c_bins = mod(g_info.c_bins-1+bin_step,g_info.nbins)+1;
g_info_r.ex_bins = mod(g_info.ex_bins-1+bin_step,g_info.nbins)+1;
end

