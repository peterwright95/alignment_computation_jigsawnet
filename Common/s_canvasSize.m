function [ out ] = s_canvasSize( orig, sz, pn, options)
%parse input
sz=double(sz);
p=sz(1); q=sz(2); m=size(orig,1); n=size(orig,2);
out = orig;
if (p~=m || q~=n)
    orig_class = class(orig);
    orig = double(orig);
    if ~exist('options','var')
        options.loc = 'center'; %expand in all directions
    end
%     if ~exist('pn','var')
%         pn = orig(1,1,:);
%     end
    pn=get_background_color(orig);
    pval=-300;
    pq=p*q;

    if (strcmp(options.loc, 'center'))
    orig_pad = padarray(orig, [floor((p-m)/2) floor((q-n)/2)], pval,'post');
    orig_pad = padarray(orig_pad, [ceil((p-m)/2) ceil((q-n)/2)], pval,'pre');
    elseif(strcmp(options.loc, 'first'))
    orig_pad = padarray(orig, [(p-m) (q-n)], pval,'post');
    end
    i = find(orig_pad(:,:,1)==-300);
    orig_pad(i) = pn(1);
    if (size(orig,3)== 3)
        orig_pad(i+pq) = pn(2);
        orig_pad(i+2*pq) = pn(3);
    end
    out=cast(orig_pad,orig_class);
end

% figure(1);imshow(lab2rgb(double(out))); title('output image');
% figure(2);imshow(lab2rgb(double(orig))); title('input image');
end


