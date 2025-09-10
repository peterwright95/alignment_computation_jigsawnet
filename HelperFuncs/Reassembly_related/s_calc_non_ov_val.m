function [ midval ] = s_calc_non_ov_val(fr_info1,fr_info2)
plt_sz=16;
alpha = 0.5;

[plt1_c, ~]=extractPalette(permute(fr_info1.c_vals_rgb,[1,3,2]),plt_sz,permute(fr_info1.c_vals,[1,3,2]));
[plt1_ex, ~ ]=extractPalette(permute(fr_info1.ex_only_vals_rgb,[1,3,2]),plt_sz,permute(fr_info1.ex_only_vals,[1,3,2]));
[plt2_c, ~]=extractPalette(permute(fr_info2.c_vals_rgb,[1,3,2]),plt_sz,permute(fr_info2.c_vals,[1,3,2]));
[plt2_ex, ~ ]=extractPalette(permute(fr_info2.ex_only_vals_rgb,[1,3,2]),plt_sz,permute(fr_info2.ex_only_vals,[1,3,2]));

dists1 = s_pdist(plt1_c,plt2_ex,'pairwise');
dists2 = s_pdist(plt2_c,plt1_ex,'pairwise');

srt_dists1 = sort(dists1(:),'ascend');
srt_dists2 = sort(dists2(:),'ascend');

midval(1) = srt_dists1(round(alpha*numel(srt_dists1)));%median(dists1(:));%
midval(2) = srt_dists2(round(alpha*numel(srt_dists2)));%median(dists2(:));%
end
