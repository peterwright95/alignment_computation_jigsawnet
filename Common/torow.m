function [ B ] = torow( A )
if ~isvector(A)
    error('Not a vector, size=(%d,%d)',size(A,1),size(A,2));
end
if iscolumn(A)
    B=A';
else
    B=A;
end
end

