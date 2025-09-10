function [ pics_info ] = s_pics_info( pics )
pattern='(?<series>[_a-zA-Z0-9]+)_p(?<part>\d+)(?<ver>_[a-zA-Z0-9]+)*\.(jpg|png)';
pics_info=regexp(pics, pattern, 'names');
end

