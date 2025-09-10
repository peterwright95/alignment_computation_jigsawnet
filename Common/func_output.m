function [ out ] = func_output( n,func,varargin )
%FUNC_OUTPUT get the n_th output variable from func
if n==2
    [~,out]=func(varargin{:});
elseif n==3
    [~,~,out]=func(varargin{:});
elseif n==4
    [~,~,~,out]=func(varargin{:});
elseif n==5
    [~,~,~,~,out]=func(varargin{:});
elseif n==6
    [~,~,~,~,~,out]=func(varargin{:});
elseif n==7
    [~,~,~,~,~,~,out]=func(varargin{:});
else
    error('not implemented for n==%d',n);
end
end

