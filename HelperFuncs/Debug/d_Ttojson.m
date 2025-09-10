function [ str ] = d_Ttojson( T,pics )
pics_info = s_pics_info(pics);
str=sprintf('{\n\t"options":\n\t{\n\t\t"series" : "%s"\n\t},\n\t"fragments":\n\t[\n\t\t',pics_info{1}.series);
for i = 1:numel(pics)
    str=[str sprintf('{\n\t\t\t"name": "%s",\n\t\t\t"T": [%0.3f, %0.3f, %0.3f]\n\t\t},',pics{i}, T(i,1),T(i,2),T(i,3))];
end
str=str(1:end-1); %remove last comma
str=[str sprintf('\n\t]\n}')];
end

