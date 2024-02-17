function [status, path] = savetiff(file, mov, ch, metadata)
%SAVETIFF This function save regmov as selected file format,
% support tiff, mat(ndmatrix)
% input:
%   - file: 1-by-1 string, the exported file name
%   - mov: 1-by-1 regmov object
%   - ch: 1-by-1 string, indicate the channel, could be "r", "g" or "b"
%   - metadata: 1-by-1 struct, with field {ftype, operation}
% output:
%   - status: the return status, SAVE_SUCCESS = 0; SAVE_FAILED = -1;
%   - path: the file saving folder
%
%   see also: bfsave

arguments
    
end
end

