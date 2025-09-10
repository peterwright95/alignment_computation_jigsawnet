function [data_mean,N,M,rdata_i,cdata_i]=s_uniform_sampling(data,res,mode)
rdata = data(:,1);
cdata = data(:,2);
rdata_low = min(rdata);
rdata_high = max(rdata);
cdata_low = min(cdata);
cdata_high = max(cdata);

dr=rdata_high-rdata_low;
dc=cdata_high-cdata_low;

norm_rdata = (rdata-rdata_low)./dr;
norm_cdata = (cdata-cdata_low)./dc;

if strcmp(mode,'old')
    N=res;
    M=res;
elseif strcmp(mode,'new')
    N=ceil(dr/res);
    M=ceil(dc/res);
else
    error('wrong mode');
end
rdata_i = floor(norm_rdata.*(N-1))+1;
cdata_i = floor(norm_cdata.*(M-1))+1;
% if any(rdata_i(:)<=0) ||any(cdata_i(:)<=0) || isempty(cdata_i) || isempty(rdata_i) ||any(rdata_i(:)>N) ||any(cdata_i(:)>M)
%     error('index 0');
% end
rdata_mean = accumarray([rdata_i,cdata_i],rdata,[N,M],@mean,NaN);
cdata_mean = accumarray([rdata_i,cdata_i],cdata,[N,M],@mean,NaN);

data_mean = round([rdata_mean(~isnan(rdata_mean)),cdata_mean(~isnan(cdata_mean))]);
end

