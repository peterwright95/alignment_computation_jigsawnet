function [nfilename]=uimwrite( varargin )
%UIMWRITE Like imwrite, but writes a unique file and does not
%overwrite the and existing file
nfilename = ufilename( varargin{2} );
varargin{2}=nfilename;
imwrite(varargin{:});
end

