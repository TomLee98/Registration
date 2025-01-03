function data = loadh5(file, opts, tspan, turbo)
%LOADH5 This function load h5 files with tspan
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
    turbo   (1,1)   logical = true %#ok<INUSA>     % use built in support
end

file = char(file);

% transform to slices span
sspan = [(tspan(1)-1)*opts.slices, tspan(2)*opts.slices-1]; % start from 0
csspan = ["/Data/Channel1/", "/Data/Channel2/"] + string((sspan(1):sspan(2))');
csspan = cellstr(csspan);

% data shape as ["X","Y","Z","T","C"]
data = zeros([opts.height, opts.width, opts.slices, diff(tspan)+1, 2], opts.dataType);

file_id = H5F.open(file);

% read channe1
for cid = 1:2
    % read slices
    for n = 1:size(csspan, 1)
        dataset_id = H5D.open(file_id, csspan{n,cid});
        zid = mod(n-1, opts.slices)+1;
        tid = ceil(n/opts.slices);
        %    X  Y   Z    T    C
        data(:, :, zid, tid, cid) = ...
            H5D.read(dataset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
        H5D.close(dataset_id);
    end
end

H5F.close(file_id);

end

