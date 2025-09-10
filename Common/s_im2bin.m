function [ bm ] = s_im2bin( rgb, lab,mode)

if strcmp(mode,'palette')
    palette_size=16;
    sigma=1.0;
    thresh=0;
    
    if ~exist('lab','var')
        lab = rgb2lab(rgb);
    end
    smooth_rgb=imgaussfilt(rgb,sigma);
    smooth_lab=imgaussfilt(lab,sigma);
    [labplette, map]=extractPalette(smooth_rgb,palette_size,smooth_lab);
    bgc = get_background_color(reshape(labplette(map(:), :),size(map,1), size(map,2),3));
    dists=pdist2(bgc,labplette);
    bg_i = find(dists<=thresh);
    mask=true(size(map));
    for i=1:numel(bg_i)
        mask=mask & bg_i(i)~=map;
    end
    mask_filled = imfill(mask, 'holes');
    mask_labels = bwlabel(mask_filled);
    lbls = unique(mask_labels);
    no_bg_lbls = lbls(lbls~=mask_labels(1));
    [~,no_bg_lbls_i] = max(arrayfun(@(lbl) sum(mask_labels(:)==lbl),no_bg_lbls));
    fr_id = no_bg_lbls(no_bg_lbls_i);
    % center = round(size(mask_filled)/2);
    % fr_id = mask_labels(center(1),center(2));
    bm = (mask_labels == fr_id);
else
    windowWidth = 11;
    polynomialOrder = 2;
    thresholdValue = 0.05;
    tol = 0.08;
    
    grayImage = lab(:,:,1)/100;%rgb(:, :, 2); % Take green channel.
    grayImage = imadjust(grayImage,stretchlim(grayImage,tol));
    % Threshold the image
    bgc = get_background_color(grayImage);
    dist_from_bg = abs(grayImage-bgc);
    binaryImage = dist_from_bg > thresholdValue;
    % Get rid of holes in the blobs.
    binaryImage = imfill(binaryImage, 'holes');
    % Extract the largest area using our custom function ExtractNLargestBlobs().
    numberToExtract = 1;
    biggestBlob = ExtractNLargestBlobs(binaryImage, numberToExtract);
    % Now get the boundaries.
    boundaries = bwboundaries(biggestBlob);
    firstBoundary = boundaries{1};
    % Get the x and y coordinates.
    x = firstBoundary(:, 2);
    y = firstBoundary(:, 1);

    % Now smooth with a Savitzky-Golay sliding polynomial filter
    smoothX = sgolayfilt(x, polynomialOrder, windowWidth);
    smoothY = sgolayfilt(y, polynomialOrder, windowWidth);
    bm=poly2mask(smoothX,smoothY,size(rgb,1),size(rgb,2));
    [bmlbls,bmlbls_n]=bwlabel(bm,4);
    if bmlbls_n>1
        bmcomp_sizes=accumarray(bmlbls(bmlbls~=0),1);
        [~,maxlbls_i] = max(bmcomp_sizes);
%         fprintf('[im2bin]-%d connected components of sizes %s detected, using component %d\n',bmlbls_n,mat2str(bmcomp_sizes),maxlbls_i);
        bm = bmlbls==maxlbls_i;
    end
end
end

%==============================================================================================
% Function to return the specified number of largest or smallest blobs in a binary image.
% If numberToExtract > 0 it returns the numberToExtract largest blobs.
% If numberToExtract < 0 it returns the numberToExtract smallest blobs.
% Example: return a binary image with only the largest blob:
%   binaryImage = ExtractNLargestBlobs(binaryImage, 1);
% Example: return a binary image with the 3 smallest blobs:
%   binaryImage = ExtractNLargestBlobs(binaryImage, -3);
function binaryImage = ExtractNLargestBlobs(binaryImage, numberToExtract)
try
	% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
	[labeledImage, numberOfBlobs] = bwlabel(binaryImage);
	blobMeasurements = regionprops(labeledImage, 'area');
	% Get all the areas
	allAreas = [blobMeasurements.Area];
	if numberToExtract > length(allAreas);
		% Limit the number they can get to the number that are there/available.
		numberToExtract = length(allAreas);
	end
	if numberToExtract > 0
		% For positive numbers, sort in order of largest to smallest.
		% Sort them.
		[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
	elseif numberToExtract < 0
		% For negative numbers, sort in order of smallest to largest.
		% Sort them.
		[sortedAreas, sortIndexes] = sort(allAreas, 'ascend');
		% Need to negate numberToExtract so we can use it in sortIndexes later.
		numberToExtract = -numberToExtract;
	else
		% numberToExtract = 0.  Shouldn't happen.  Return no blobs.
		binaryImage = false(size(binaryImage));
		return;
	end
	% Extract the "numberToExtract" largest blob(a)s using ismember().
	biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
	% Convert from integer labeled image into binary (logical) image.
	binaryImage = biggestBlob > 0;
catch ME
	errorMessage = sprintf('Error in function ExtractNLargestBlobs().\n\nError Message:\n%s', ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
end
% % img   : the image to binarize
% % method: the gradient method used
%     %----Contants-------
%     DEBUG_VERB = 0;
%     N_FIG = 200;
%     %-------------------
%     if ~exist('method','var')
%         method='sobel';
%     end
%     if ~exist('cs','var')
%         cs='rgb';
%     end
%     [I,~] = imgradient(s_im2gray(img, cs), method);
%     I = I / max(I(:));
%     if (DEBUG_VERB >=3)
%         figure(N_FIG), imshow(I), title('[Im2Bin]: Gradient Map');
%         N_FIG = N_FIG + 1;
%     end
%     BW = imbinarize(I, 'adaptive'); %available only in 2016 version
% %     if contains(version,'2017')
% %         BW = imbinarize(I, 'adaptive'); %available only in 2016 version
% %     else
% %         BW=adaptivethreshold(I,15,0.01,0);
% %         BW=~imfill(~BW,[1 1]);
% %     end    
%     if (DEBUG_VERB >=2)
%         figure(N_FIG), imshow(~BW), title('[Im2Bin]: Adaptive binarization');
%         N_FIG = N_FIG + 1;
%     end
%     BW = imfill(BW, 'holes');
% 
%     CC = bwconncomp(BW, 4);
%     S = regionprops(CC, 'Area');
%     L = labelmatrix(CC);
%     [~,maxi] = max([S.Area]);
%     bm = ismember(L, maxi);
%     if (DEBUG_VERB >=1)
%         figure(N_FIG), imshow(bm), title('[Im2Bin]: Final Binary Image');
%         N_FIG = N_FIG + 1;
%     end
%     
%     
% %     [I,~] = imgradient(s_im2gray(img, cs), method);
% %     %[I,~,~] = obj.grad(img);
% %     maxI = max(I(:));
% %     if (maxI ~= 0)
% %         I = I / maxI;
% %         %I(I>0) = 1;
% %         if (DEBUG_VERB >=3)
% %             figure(N_FIG), imshow(I), title('[Im2Bin]: Gradient Map');
% %             N_FIG = N_FIG + 1;
% %         end
% %         if (any(isnan(I(:))))
% %             error('gradient contains nan');
% %         end
% %         complete = imfill(I, 'holes');
% %         if (DEBUG_VERB >=2)
% %             figure(N_FIG), imshow(complete), title('[Im2Bin]: No holes');
% %             N_FIG = N_FIG + 1;
% %         end
% % %         [his,bin] = imhist(complete);
% % %         nhis = his / numel(complete);
% % %         chis = cumsum(nhis);
% % %         chis = chis - nhis(1); % first normalized count  belongs to the 0 value bin which we dont need
% % %         thresh = bin(find(chis>0.1,1)); % look for the first non zero intensity which acumulates to 2% small values
% % %         complete(complete < thresh) = 0;
% %         complete(complete < 0.1) = 0;
% %         bm = logical(complete);
% %         if (DEBUG_VERB >=1)
% %             figure(N_FIG), imshow(bm), title('[Im2Bin]: Final Binary Image');
% %             N_FIG = N_FIG + 1;
% %         end
% %     else
% %         bm = I;
% %     end


