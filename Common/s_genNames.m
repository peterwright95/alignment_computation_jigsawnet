function [ pics ] = s_genNames( series_name, parts_num, postfix, ext)
if ~exist('postfix','var')
    postfix='';
end
if ~exist('ext','var')
    ext='.jpg';
end
pics=arrayfun(@(x) {[series_name, '_p',num2str(x),postfix,ext]},parts_num);
end

