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
tsspan = sspan + 1;     % start from 1
tsspan = ceil((tsspan(1):tsspan(2))./opts.slices);
csspan = ["/Data/Channel1/", "/Data/Channel2/"] + string((sspan(1):sspan(2))');

% data shape as ["X","Y","Z","T","C"]
data = zeros([opts.height, opts.width, opts.slices, diff(tspan)+1, 2], opts.dataType);

% read channe1
for c = 1:2
    % read slices
    for n = 1:size(csspan, 1)
        %    X  Y        Z          T       C
        data(:, :, sspan(1)+n, tsspan(n), c) = h5read(file, csspan(n,c));
    end
end

end

