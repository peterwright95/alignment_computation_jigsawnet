function d=s_tform_dist(t1,t2)
%Calculates the transformation distance between the set t1 of
%transformations to the set t2 of transformations.
%t1 is of size m x 3, t2 is of size n x 3, where m is the number of
%transformations in the set t1 and n is the number of transformations in t2
%

if iscolumn(t1)
    t1=t1';
end
if iscolumn(t2)
    t2=t2';
end
if size(t1,2)~=3 || size(t2,2) ~= 3
    error('prog:s_tform_dist',...
          'Wrong tform set size: size(t1)=[%d,%d], size(t2)=[%d,%d]\n(expected size(.,2)=3, size(t1,1)=size(t2,1))',...
          size(t1,1),size(t1,2),size(t2,1),size(t2,2));
end
if size(t1,1)==1
    t1 = repmat(t1,size(t2,1),1);
elseif size(t2,1)==1
    t2 = repmat(t2,size(t1,1),1);
end
roterr=acos(cosd(t1(:,3)-t2(:,3)))/pi;

trerr=sqrt(sum((t1(:,1:2)-t2(:,1:2)).^2,2));
trerr=trerr./(75+trerr);

d=trerr+sqrt(roterr);
end

