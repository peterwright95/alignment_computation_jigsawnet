function [dist] = s_polys_dist(poly1,poly2,extrap_intersectionPol)
DEBUG_VERB = 0;
is_restrict_to_extrapolated_region = exist('extrap_intersectionPol','var');
[V1,E1]=cell2poly(poly1);
[V2,E2]=cell2poly(poly2);
if is_restrict_to_extrapolated_region
    [V1_2_e,E1_2_e]=cell2poly(extrap_intersectionPol);
end
D=pdist2(V1,V2);
[~,minD_p1]=min(D,[],2);
[~,minD_p2]=min(D,[],1);
bb_i1=minD_p2(minD_p1);
if ~isrow(bb_i1)
    bb_i1=bb_i1';
end
bb_i2=minD_p1(minD_p2);
if ~isrow(bb_i2)
    bb_i2=bb_i2';
end

p1_bbi=find(bb_i1==1:numel(bb_i1)); % poly 1 best buddy index
p2_bbi=minD_p1(p1_bbi)'; % poly 2 best buddy index
if numel(p1_bbi) ~= numel(p2_bbi)
    error('number of best buddies should be equal. n_bb1=%d, n_bb2=%d',numel(p1_bbi),numel(p2_bbi));
end
unordered_p2_bbi=find(bb_i2==1:numel(bb_i2));
if ~all(ismember(p2_bbi,unordered_p2_bbi))
    error('poly2 best buddy is not the same as through poly1');
end
if is_restrict_to_extrapolated_region
    isV1inExtP=inpoly(V1(p1_bbi,:),V1_2_e,E1_2_e);
    isV2inExtP=inpoly(V2(p2_bbi,:),V1_2_e,E1_2_e);
    isInExtP = isV1inExtP&isV2inExtP;
else
    isInExtP = true(size(p1_bbi));
end
p1_bbi_inext = p1_bbi(isInExtP);
p2_bbi_inext = p2_bbi(isInExtP);

if ~any(isInExtP)
    p1_bbi_inext = p1_bbi;
    p2_bbi_inext = p2_bbi;
end
[isV1inP2,isV1onP2]=inpoly(V1(p1_bbi_inext,:),V2,E2);
[isV2inP1,isV2onP1]=inpoly(V2(p2_bbi_inext,:),V1,E1);
isV1strictly_inP2 = isV1inP2 & ~isV1onP2;
isV2strictly_inP1 = isV2inP1 & ~isV2onP1;
dist_bbi = sub2ind(size(D),p1_bbi_inext,p2_bbi_inext);
signed_D = D(dist_bbi);
turn_sign_cond = isV2strictly_inP1|isV1strictly_inP2;
signed_D(turn_sign_cond) = -signed_D(turn_sign_cond);
[~,idx]=max(abs(signed_D));
dist = signed_D(idx);
% dist = median(signed_D);
% dist = mean(signed_D);
% dist = min(signed_D);
% if dist>=0
%     dist = max(signed_D);
% end

if DEBUG_VERB >=1
    figure(1);
    plot(V1(:,2),V1(:,1),V2(:,2),V2(:,1)...
        ,V1(p1_bbi_inext(~isV1strictly_inP2),2),V1(p1_bbi_inext(~isV1strictly_inP2),1),'+',V2(p2_bbi_inext(~isV2strictly_inP1),2),V2(p2_bbi_inext(~isV2strictly_inP1),1),'+'...
        ,V1(p1_bbi_inext(isV1strictly_inP2),2),V1(p1_bbi_inext(isV1strictly_inP2),1),'o',V2(p2_bbi_inext(isV2strictly_inP1),2),V2(p2_bbi_inext(isV2strictly_inP1),1),'o');
    if is_restrict_to_extrapolated_region
    hold on;
    plot(V1_2_e(:,2),V1_2_e(:,1)...
        ,V1(p1_bbi(~isInExtP),2),V1(p1_bbi(~isInExtP),1),'*',V2(p2_bbi(~isInExtP),2),V2(p2_bbi(~isInExtP),1),'*');
    hold off;
    end
    title(sprintf('Distance=%0.3f',dist));
end


% [d_min1, x_d_min1, y_d_min1] = p_poly_dist(poly1(:,1)', poly1(:,2)', poly2(:,1)', poly2(:,2)');
% [min_dist1,min_dist1_i] = min(d_min1);
% [d_min2, x_d_min2, y_d_min2] = p_poly_dist(poly2(:,1)', poly2(:,2)',poly1(:,1)', poly1(:,2)');
% [min_dist2,min_dist2_i] = min(d_min2);
% 
% if min_dist1<min_dist2
%     p_in_poly1 = poly1(min_dist1_i,:);
%     p_in_poly2 = [x_d_min1(min_dist1_i),y_d_min1(min_dist1_i)];
%     dist = min_dist1;
% else
%     p_in_poly2 = poly2(min_dist2_i,:);
%     p_in_poly1 = [x_d_min2(min_dist2_i),y_d_min2(min_dist2_i)];
%     dist = min_dist2;
% end


end

function [verts,edges] = cell2poly(polyscell)
if iscell(polyscell)
    if size(polyscell,2)~= 1
        polyscell = polyscell';
    end
    isNx2=cellfun(@(x) size(x,2) == 2,polyscell);
    if any(~isNx2)
        polyscell(~isNx2) = cellfun(@(x) {x'},polyscell(~isNx2));
    end
    verts = cell2mat(polyscell);
    n_verts = cellfun(@(x) size(x,1),polyscell);
    endpoints=cumsum(n_verts);
    startpoints=endpoints-n_verts+1;
    edges = cell2mat(arrayfun(@(s,e) {[(s:e-1)',(s+1:e)';e,s]},startpoints,endpoints));
%     prevcell = cell(numel(polyscell),1);
%     prevcell(2:end) = polyscell(1:end-1);
%     edges = cell2mat(cellfun(@(x,y) {poly_edges(x)+size(y,1)},polyscell,prevcell));
else
    if size(polyscell,2) ~= 2
        polyscell = polyscell';
    end
    verts = polyscell;
    edges = poly_edges(polyscell);
end
end

function edges=poly_edges(poly)
    n_vert = size(poly,1);
    edges = [(1:n_vert-1)',(2:n_vert)';...
             n_vert, 1];
end




