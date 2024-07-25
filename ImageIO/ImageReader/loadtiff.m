function data = loadtiff(file, opts, tspan, turbo)
%LOADTIFF This function load tiff files with tspan
% Input:
%   - file: 1-by-1 string, the data full file path
%   - opts: 1-by-12 table, images properties
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - turbo: 1-by-1 logical, use external library for fast loading
% Output:
%   - info: struct with {opts, rt}, where opts is 1-by-12 metadata table,
%     rt is n-by-1 nonnegtive double time array
%
% where you can pass the time span as tspan for continous time cut

arguments
    file    (1,1)   string
    opts    (1,12)  table
    tspan   (1,2)   double  {mustBePositive, mustBeInteger}
    turbo   (1,1)   logical = false
end

file = char(file);

if turbo == true
    data = pyopen_reg(file);
else
    % transform tspan to sspan
    if opts.dimOrder(end) == "T"
        sspan = [(tspan(1)-1)*opts.slices*opts.channels+1, ...
            tspan(2)*opts.slices*opts.channels];
        data = bfopen_reg(file, sspan);
    else
        throw(MException("loadtiff:invalidDimension", ...
            "Dimension <T> is not the last."));
    end
end

end

function mov = pyopen_reg(file)
% check the python environment path
pathToReader = fileparts(mfilename('fullpath'));
pathToLTPY = [pathToReader, filesep, 'load_tiff.py'];
if count(py.sys.path, pathToLTPY) == 0
    insert(py.sys.path,int32(0), pathToLTPY);
end
fname = py.str(file);

try
    % imagej stack: 'TZCYXS'
    mov = pyrunfile("load_tiff.py", "vol", file=fname);
catch ME
    throwAsCaller(ME);
end
% convert img from ndarray to matlab value
mov = cast(mov.astype('single'), "uint16");
if ndims(mov) == 3
    mov = permute(mov, [2,3,1]);   % to (Y,X,Z)
elseif ndims(mov) == 4
    mov = permute(mov, [3,4,2,1]); % to (Y,X,Z,T)
elseif ndims(mov) == 5
    mov = permute(mov, [4,5,3,2,1]); % to (Y,X,C,Z,T)
else
    throw(MException("imload:invalidImagesStackDimension", ...
        "Image stack dimension > 5 is not supported."));
end
end