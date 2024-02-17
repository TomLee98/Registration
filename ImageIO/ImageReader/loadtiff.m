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
    tspan   (1,:)   double = []
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
    end
end

% reconstruct the image stack
tmpopts = opts;
tmpopts.frames = diff(tspan)+1;
data = imreshape(data, tmpopts);

end

function mov = pyopen_reg(file)
% check the python environment path
if count(py.sys.path,'/ImageIO/ImageReader/load_tiff.py') == 0
    insert(py.sys.path,int32(0), ...
        '/ImageIO/ImageReader/load_tiff.py');
end
fname = py.str(file);

try
    % imagej stack: 'TZCYXS'
    mov = pyrunfile("load_tiff.py", "vol", file=fname);
catch ME
    throwAsCaller(ME);
end
% convert img from ndarray to matlab value
mov = cast(mov.astype('single'), dataType);
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