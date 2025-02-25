function data = loadnd2(file, opts, tspan, turbo)
%LOADND2 This function load nd2 files with tspan
% Input:
%   - file: 1-by-1 string, the data full file path
%   - opts: 1-by-12 table, images properties
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - turbo: 1-by-1 logical, use external library for fast loading
% Output:
%   - data: n dimension MATLAB array
%
% where you can pass the time span as tspan for continuous time cut

arguments
    file    (1,1)   string
    opts    (1,12)  table
    tspan   (1,2)   double  {mustBePositive, mustBeInteger}
    turbo   (1,1)   logical = false
end

file = char(file);

if opts.dimOrder(end) == "T"
    
    if turbo == true
        sspan = [(tspan(1)-1)*opts.slices+1, ...
            tspan(2)*opts.slices];
        data = nd2open_reg(file, sspan);
    else
        sspan = [(tspan(1)-1)*opts.slices*opts.channels+1, ...
            tspan(2)*opts.slices*opts.channels];
        data = bfopen_reg(file, sspan);
    end
else
    throw(MException("loadnd2:invalidDimension", ...
        "Dimension <T> is not the last."));
end

end

function mov = nd2open_reg(file, sspan)
% use nd2 library load total volumes
[FilePointer, ImagePointer, ImageReadOut] = ND2Open(file);
% mov is a c-by-1 cell, c for channels number
mov = ND2Read(FilePointer, ImagePointer, ImageReadOut, sspan(1):sspan(2));
ND2Close(FilePointer)
clear("FilePointer","ImagePointer","ImageReadOut");

for n = 2:numel(mov)
    mov{1} = cat(4, mov{1}, mov{n});
    mov{n} = [];    % free memory
end
mov = mov{1};
% permute as (Y,X,C,Z(T)): nd2 file fixed dimension order
mov = permute(mov, [1,2,4,3]);
end