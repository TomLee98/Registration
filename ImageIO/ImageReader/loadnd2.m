function data = loadnd2(file, opts, tspan, turbo)
%LOADND2 This function load nd2 files with tspan
% Input:
%   - file: 1-by-1 string, the data full file path
%   - opts: 1-by-12 table, images properties
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - turbo: 1-by-1 logical, use external library for fast loading
% Output:
%   - data: ndimension matlab array
%
% where you can pass the time span as tspan for continous time cut

arguments
    file    (1,1)   string
    opts    (1,12)  table
    tspan   (1,:)   double  {mustBeVector}
    turbo   (1,1)   logical = false
end

file = char(file);

if turbo == true
    data = nd2open_reg(file, tspan);
else
    if opts.dimOrder(end) == "T"
        sspan = [(tspan(1)-1)*opts.slices*opts.channels+1, ...
            tspan(2)*opts.slices*opts.channels];
        data = bfopen_reg(file, sspan);
    else
        throw(MException("loadnd2:invalidDimension", ...
            "Dimension <T> is not the last."));
    end
end

% reconstruct the image stack
tmpopts = opts;
tmpopts.frames = diff(tspan)+1;
data = imreshape(data, tmpopts);

end

function mov = nd2open_reg(file, tspan)
% use nd2 library load total volumes
[FilePointer, ImagePointer, ImageReadOut] = ND2Open(file);
% mov is a c-by-1 cell, c for channels number
mov = ND2Read(FilePointer, ImagePointer, ImageReadOut, tspan(1):tspan(2));
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

function mov = imreshape(mov, opts)
bindingArray = ["Y",    "X",    "C",    "Z",    "T";...
    "width","height","channels","slices","frames"];
[~,Loc] = ismember(opts.dimOrder,bindingArray(1,:));
mov = reshape(mov,...
    opts.(bindingArray(2,Loc(1))),...
    opts.(bindingArray(2,Loc(2))),...
    opts.(bindingArray(2,Loc(3))),...
    opts.(bindingArray(2,Loc(4))),...
    opts.(bindingArray(2,Loc(5))));
end