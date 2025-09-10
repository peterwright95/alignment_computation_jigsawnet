function [ invtforms ] = s_invtform( tforms )
if size(tforms,2) ~= 3
    tforms = tforms';
    if size(tforms,2) ~= 3
        error('tforms must be of dimension 3');
    end
end
invtforms = -tforms;
R1 = [cosd(invtforms(:,3)),-sind(invtforms(:,3))];
R2 = [-R1(:,2),R1(:,1)];
invtforms(:,1:2) = round([sum(R1.*invtforms(:,1:2),2),sum(R2.*invtforms(:,1:2),2)]);
end

