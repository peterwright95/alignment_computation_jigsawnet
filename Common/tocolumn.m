function [ B ] = tocolumn( A )
if isempty(A)
    B=A;
else
    if ~isvector(A)
        error('Not a vector, size=(%d,%d)',size(A,1),size(A,2));
    end
    if isrow(A)
        B=A';
    else
        B=A;
    end
end
end

