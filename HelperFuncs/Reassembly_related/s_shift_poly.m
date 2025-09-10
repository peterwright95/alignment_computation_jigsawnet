function spoly=s_shift_poly(poly,t)
R=rot2D(t(3));
% spoly = round(R'*poly + repmat(flip(t(1:2)'),1,size(poly,2)));
spoly = round(bsxfun(@plus,R'*poly,flip(t(1:2)')));
end
