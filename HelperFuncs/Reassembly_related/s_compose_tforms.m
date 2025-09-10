function [ outtforms] = s_compose_tforms( tforms1,tforms2 )
%S_COMPOSE_TFORMS - create the transforms composition of T1(T2)
%for each pair of transforms T1=tforms1(i,:), and T2=tforms(i,:)
%
if ~isequal(size(tforms1),size(tforms2))
    error('input must be the same size');
end
if size(tforms1,2) ~= 3
    tforms1 = tforms1';
    tforms2 = tforms2';
    if size(tforms1,2) ~= 3
        error('tforms must be of dimension 3');
    end
end

outtforms = zeros(size(tforms1));
outtforms(:,3) = tforms1(:,3) + tforms2(:,3);
R1 = [cosd(tforms1(:,3)),-sind(tforms1(:,3))];
R2 = [-R1(:,2),R1(:,1)];
outtforms(:,1:2) = round([sum(R1.*tforms2(:,1:2),2),sum(R2.*tforms2(:,1:2),2)] + tforms1(:,1:2));
end

