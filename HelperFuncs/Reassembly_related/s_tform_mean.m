function [ mean_tform ] = s_tform_mean( tforms )
n_tforms=size(tforms,1);
if n_tforms > 1
    if size(tforms,2)~=3
        error('size 2nd dimension of tforms should equal to 3, instead was %d',size(tforms,2));
    end
    mean_tform = zeros(1,3);
    R=zeros(3);
    for i=1:n_tforms
        R=R+rotz(tforms(i,3));
    end
    [U,~,V]=svd(R);
    V=V';
    S=diag([1,1,det(U)*det(V)]);
    meanR=U*S*V;
    mean_r=rotm2eul(meanR);
    mean_tform(3)=radtodeg(mean_r(1));
    mean_tform(1:2)=mean(tforms(:,1:2),1);
else
    mean_tform=tforms;
end
end

