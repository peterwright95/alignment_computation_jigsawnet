function [ cspace,cspace_ex] = s_calc_allcspace( fr_input,pairs)
    P=cellfun(@(x) x.poly(:,1:end-1),fr_input,'UniformOutput',false);
    P_ex=cellfun(@(x) x.poly_ex(:,1:end-1),fr_input,'UniformOutput',false);
    tic;
    cspace = calc_cspace(P,{},pairs',0:359);
    fprintf('All cspace calc time(mex): %0.4f sec\n',toc);    
    tic;
    cspace_ex = calc_cspace(P_ex,P,pairs',0:359);
    fprintf('All cspace (ex) calc time(mex): %0.4f sec\n',toc);
end