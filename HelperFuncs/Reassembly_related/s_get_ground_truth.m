function [ absT,relT ] = s_get_ground_truth( seriesname,pairs )
injson_filename = 'groundtruth.json';
fileID = fopen(injson_filename,'r');
st=jsondecode(fscanf(fileID,'%s'));
fclose(fileID);
absT = [];
relT = [];
for i=1:numel(st)
    options_st = st(i).options;
    if (strcmp(options_st.series,seriesname))
        frags_st = st(i).fragments;
        pics_info = s_pics_info({frags_st.name});
        parts = cellfun(@(x) str2double(x.part),pics_info);
        absT = zeros(max(parts),3);
        absT(parts,:) = cell2mat({frags_st.T})';
        break;
    end
end
if ~isempty(absT) && exist('pairs','var')
    tform_p1 = absT(pairs(:,1),:);
    tform_p2 = absT(pairs(:,2),:);
    relT = s_compose_tforms(s_invtform(tform_p1),tform_p2);
end
end

