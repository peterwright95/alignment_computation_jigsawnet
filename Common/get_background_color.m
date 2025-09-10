function bg=get_background_color(m)
corners=reshape([m(1,1,:);m(1,end,:);m(end,1,:);m(end,end,:)],4,size(m,3));
[uc,~,ui]=unique(corners,'rows');
counter = accumarray(ui,1,[size(uc,1),1]);
[~,maxi]=max(counter);
bg=uc(maxi,:);
end
