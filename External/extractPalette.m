function [labplette, pxl2palette_map, palette, palette_im] = extractPalette(rgb, palette_size, lab)
  if nargin < 1
      display('please input rgb image, default palette size is 5')
  else
      if nargin < 2
         display('default palette size is 5')
         palette_size = 5;
      elseif nargin < 3
          lab = rgb2lab(rgb);
      end
      if nargout >= 3
          [labplette, pxl2palette_map,palette] = calc_palette(rgb,lab,palette_size);
          if (nargout >= 4)
              PAL_IM_SIZE = 100;

              palette_im=reshape(palette, 1,palette_size,3);
              palette_im=repmat(palette_im,PAL_IM_SIZE*PAL_IM_SIZE,1,1);
              palette_im=uint8(reshape(palette_im,PAL_IM_SIZE,PAL_IM_SIZE*palette_size,3));
          end
      else
          [labplette, pxl2palette_map] = calc_palette(rgb,lab,palette_size);
      end
  end
  
function [ labpalette, pxl2palette_map, palette ] = calc_palette( rgb, lab, palette_size )
    [height, width, channel] = size(rgb);
    dim = double(rgb);
    assert(channel == 3);
    ngrid = 16;
    step_size = 255 / (ngrid - 1);
    sample_bin = zeros(ngrid*ngrid*ngrid, 4);
    
    bin = round(dim / step_size);
    bin_num_full = bin(:,:,1) * (ngrid*ngrid) + bin(:,:,2)*ngrid + bin(:,:,3) + 1;
    bin_num = bin_num_full(:);
    maxidx = max(bin_num);
    for d=1:3
        sample_bin(1:maxidx,d) = accumarray(bin_num, reshape(lab(:,:,d),height*width,1));
    end
    sample_bin(1:maxidx,4) = accumarray(bin_num, ones(height*width,1));
    nz_bin_i = find(sample_bin(:,4)~=0);
    sample_bin_not_empty = sample_bin(nz_bin_i,:);
    D = bsxfun(@rdivide, sample_bin_not_empty(:,1:3), sample_bin_not_empty(:,4));
    sample_cnt = sample_bin_not_empty(:,4);
    
    %Initialize palette
    pickcnt = sample_cnt;
    palette = zeros(palette_size+1, 3);
    for ii = 1 : palette_size
        [~, idx] = max(pickcnt);
        palette(ii, :) = D(idx, :);
        dis = sum(bsxfun(@minus, D(idx,:), D).^2, 2) / (80*80);
        pickcnt = pickcnt .* (1-exp(-dis));
    end
    
    % add black
        
    szD = size(D,1);
    szPlt = palette_size+1;
    for ii = 1: 20
        sumD = zeros(palette_size+1, 3);
                
        Di = repmat((1:szD)', 1, szPlt);
        Plti = repmat((1:szPlt), szD,1);
        a=D(Di(:),:);
        b=palette(Plti(:),:);
        [~,minI] = min(reshape(sum((a-b).^2,2),szD,szPlt),[],2); %get the closest palette color for each bin representative
        cnt = accumarray(minI, sample_cnt);
        for d=1:3
            sumD(:,d) = accumarray([minI;(1:szPlt)'], [sample_cnt.*D(:,d);zeros(szPlt,1)]);
        end
        cnt_i = find(cnt>0);
        cnt_i(cnt_i==szPlt) = [];
        palette(cnt_i,:) = bsxfun(@rdivide, sumD(cnt_i,:), cnt(cnt_i));
    end
    labpalette = palette(1:palette_size,:);
    if (nargout>=3)
        palette = lab2rgb(labpalette,'OutputType','uint8');
    end
    % pixel to palette mapping
    bin2palette_map = zeros(ngrid*ngrid*ngrid, 1);
    bin2palette_map(nz_bin_i) = minI;
    pxl2palette_map = bin2palette_map(bin_num_full);
    
    blk = repmat([0 0 0],palette_size,1);
    [~,closest2blk] = min(sum((labpalette-blk).^2,2));
    pxl2palette_map(pxl2palette_map==szPlt) = closest2blk;
end
end