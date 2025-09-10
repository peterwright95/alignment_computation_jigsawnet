function [ cycdiff ] = cyclic_diff(a,b)
%CYCLIC_DIFF Summary of this function goes here
%   Detailed explanation goes here
if ~exist('b','var')
    abdiff = diff(a,1,2);
else
    if ~isequal(size(a),size(b))
        error('sizeof a and b must be equal');
    end
    if any(b<0|b>1)
        error('b values must be in [0,1] interval');
    end
    abdiff = b-a;
end
if any(a<0|a>1)
    error('a values must be in [0,1] interval');
end

abdiff = mod(abdiff,1+eps);
cycdiff = abdiff;
% cycdiff = min(abdiff,1-abdiff);
end

