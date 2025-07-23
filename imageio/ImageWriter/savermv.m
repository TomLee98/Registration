function [status, path] = savermv(file, mov, ch, tspan, metadata, turbo)
%SAVEMAT This function save regmov as mat file format
% input:
%   - file: 1-by-1 string, the exported file name
%   - mov: 1-by-1 regmov object
%   - ch: 1-by-c string, indicate the channel, could be "r", "g" or "b"
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - metadata: 1-by-1 struct, with field {ftype, operation}
%   - %   - turbo: 1-by-1 logical, use external library for fast saving
% output:
%   - status: the return status, SAVE_SUCCESS = 0; SAVE_FAILED = -1;
%   - path: the file saving folder
%
%   see also: regmov

arguments
    file        (1,1)   string
    mov         (1,1)   regmov
    ch          (1,:)   string  {mustBeMember(ch, ["r","g","b"])}
    tspan       (1,2)   double  {mustBePositive, mustBeInteger} %#ok<INUSA> fully save
    metadata    (1,1)   struct
    turbo       (1,1)   logical = true %#ok<INUSA>         % no matter 'turbo' is
end

% memory size validation
mem = GetAvailableMemory();
if mem < (mov.Bytes.mem + mov.Bytes.map)
    warning("No enough memory for tiff saving.");
    status = -1;
    return;
end

[path, ~, ~] = fileparts(file);

try
    Time = mov.Time;

    MetaData = mov.MetaData;

    % extract data from regmov to memory
    chs = MetaData.cOrder;
    Data = mov.Movie(:,:,ch==chs,:,:);   % 5D array
    MetaData.cOrder = ch;       % modify the output color channel
    MetaData.channels = numel(ch);

    Supplementary = metadata;

    % no compression for fast saving
    save(file, "Data", "MetaData", "Supplementary", "Time", "-v7.3", ...
        "-nocompression");
    
    status = 0;
catch ME
    rethrow(ME);
end


end

