function A = s_imshift(im, shift, outputview)
if (nargin<3)
    outputview='same';
end
shiftRows = round(shift(1)); shiftCols = round(shift(2)); 
if (strcmp(outputview,'same'))
    A = repmat(im(end,1,:),size(im,1),size(im,2));

    if shiftRows >= 0 && shiftCols >= 0
        A(1+shiftRows:end,1+shiftCols:end,:) = im(1:end-shiftRows,1:end-shiftCols,:);
    elseif shiftRows >= 0 && shiftCols < 0
        A(1+shiftRows:end,1:end+shiftCols,:) = im(1:end-shiftRows,1-shiftCols:end,:);
    elseif shiftRows < 0 && shiftCols >= 0
        A(1:end+shiftRows,1+shiftCols:end,:) = im(1-shiftRows:end,1:end-shiftCols,:);
    else
        A(1:end+shiftRows,1:end+shiftCols,:) = im(1-shiftRows:end,1-shiftCols:end,:);
    end
else
    M=size(im,2);
    if (shiftRows<0)
        A = [im;repmat(im(end,1,:),-shiftRows,M)];
    else
        A = [repmat(im(end,1,:),shiftRows,M); im];
    end
    N = size(A,1);
    if (shiftCols<0)
        A = [A,repmat(im(end,1,:),N,-shiftCols)];
    else
        A = [repmat(im(end,1,:),N,shiftCols),A];
    end
end
end
