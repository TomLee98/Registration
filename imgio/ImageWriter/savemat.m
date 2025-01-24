function [status, path] = savemat(file, mov, ch, tspan, metadata, turbo)
%SAVEMAT This function save regmov as mat file format
% input:
%   - file: 1-by-1 string, the exported file name
%   - mov: 1-by-1 regmov object
%   - ch: 1-by-1 string, indicate the channel, could be "r", "g" or "b"
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - metadata: 1-by-1 struct, with field {ftype, operation}
%   - %   - turbo: 1-by-1 logical, use external library for fast saving
% output:
%   - status: the return status, SAVE_SUCCESS = 0; SAVE_FAILED = -1;
%   - path: the file saving folder
%
%   see also: bfsave

status = 0;
path = "";
end

