function [ ufilename ] = ufilename( filename )
%UFILENAME find a unique filename. also create dir if does not exist

[dir,name,ext] = fileparts(filename);
if ext(1) ~= '.'
    error('wrong extension extraction from filename');
end
ext=ext(2:end);
if (exist(dir, 'dir') == 0)
    mkdir(dir);
end
filesCount = 0;
ufilename = fullfile(dir, [name '.' ext]);
while(exist(ufilename, 'file') == 2)
    filesCount = filesCount + 1;
    ufilename = fullfile(dir, [name '_' num2str(filesCount) '.' ext]);
end
end

