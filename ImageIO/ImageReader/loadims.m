function data = loadims(file, opts, tspan, turbo)
%LOADIMS This function load ims files with tspan
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
    % generate HDF5 reader
    GID = H5G.open(H5F.open(file), '/DataSet');
    dataobj = DatasetReader(GID);

    % load the image data, X,Y,Z,C,T
    data = dataobj.GetData(tspan);
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

end
