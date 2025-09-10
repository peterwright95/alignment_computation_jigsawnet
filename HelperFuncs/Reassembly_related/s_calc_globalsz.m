function globalsz=s_calc_globalsz(fr_input,pairs)
    n_pairs = size(pairs,1);
    globalsz=zeros(n_pairs,1);
    for ci=1:n_pairs
        p1=pairs(ci,1);
        p2=pairs(ci,2);
        in_p1=fr_input{p1};
        in_p2=fr_input{p2};
        diagf_l     = norm(in_p1.origsz_ex);
        diagnf_l    = norm(in_p2.origsz_ex);
        globalsz(ci) = max(round(diagf_l + 2*diagnf_l),round(2*diagf_l + diagnf_l));
    end
end

