function [ Si,Mi ] = s_reorder_mask( mask )
%S_REORDER_MASK find the cyclic order of the mask contour
%   input: 
%       mask - NxM logical matrix of a fragment
%   output:
%       Si - the cyclic order of the contour points indices.
%            Namely, Si(1)=index of first mask contour point in order.
%            Si(2)=index of second mask contour point in order...
%       Mi - the cyclic order location for each contour point.
%            Namely, Mi(1)=the cyclic order location of contour point 1
%            Mi(2)=the cyclic order location of contour point 2...
%
%
%
M=s_mask_contour(mask);
n=max(size(M));
M=s_canvasSize(M,[n,n]);
DblM = false(n*2,n);
DblM(1:2:n*2,:) = M;
% bb=get_bounding_box(M);
% DblM(2*bb.rmin-1,:)=false; % disable first row to avoid choosing it as cutting points.
% DblM(2*bb.rmax-1,:)=false; % disable last row to avoid choosing it as cutting points.

lbls = bwlabel(DblM);
[~,slbls_i] = sort(accumarray(lbls(lbls~=0),1),'ascend');
for i=1:numel(slbls_i)
    lbls_i=slbls_i(i);
    [cutp_r,cutp_c] = find(lbls==lbls_i);
    cutp_r = (cutp_r+1)./2;
    if func_output(2,@bwlabel,M(cutp_r,:))>1
        cutp_i = sub2ind([n,n],cutp_r,cutp_c);
        cutp_c_min = min(cutp_c);
        cutp_c_max = max(cutp_c);
        if any(M(cutp_r(1)-1,cutp_c_min-1:cutp_c_max+1)) && ...
           any(M(cutp_r(1)+1,cutp_c_min-1:cutp_c_max+1))
           break; 
        end
    end
end
curp_mask = false(n);
curp_mask(cutp_r(1)-1:cutp_r(1)+1,cutp_c_min-1:cutp_c_max+1)=[true(1,cutp_c_max-cutp_c_min+3);...
                       false(1,cutp_c_max-cutp_c_min+3);...
                       true(1,cutp_c_max-cutp_c_min+3)];
[start_points(1,:),start_points(2,:)]=find(curp_mask&M);
sp=start_points(:,1);
W = ones(n);% use constant metric
L = -Inf(n);% create a mask to restrict propagation
L(M) = +Inf;
L(cutp_i)=-Inf;
options.constraint_map = L;
D = perform_fast_marching(W, sp, options);% do the FM computation

[mask_rc,mask_i]=mask_indices(M);
inf_pixels_i = find(isinf(D(M)));
infmax = @(x) max(x(~isinf(x)));
for i=1:numel(inf_pixels_i)
    inf_rc = mask_rc(inf_pixels_i(i),:);
    N_Di=D((inf_rc(1)-1):(inf_rc(1)+1),(inf_rc(2)-1):(inf_rc(2)+1));
    if all(isinf(N_Di(:)))
        error('pixel %s has no close non inf pixels after fase marching',mat2str(inf_rc));
    else
        D(mask_i(inf_pixels_i(i)))=infmax(N_Di);
    end    
end

[~,Si]=sort(D(M));
Mi(Si)=(1:numel(Si))./numel(Si);

% figure(1);subplot(1,2,1);imshow(convert_distance_color(D)); subplot(1,2,2);[tmpM(:,1),tmpM(:,2)]=find(M);plot(tmpM(Si,2),-tmpM(Si,1));axis equal;
end

