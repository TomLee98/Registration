function data = loadrmv(file, opts, tspan, turbo)
%LOADMAT This function load mat files with tspan
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
    opts    (1,12)  table %#ok<INUSA>           
    tspan   (1,2)   double  {mustBePositive, mustBeInteger} %#ok<INUSA> fully loading
    turbo   (1,1)   logical = true             %#ok<INUSA> % no matter 'turbo' is
end

try
    load(file, "-mat", "Data");
    data = Data;
catch ME
    rethrow(ME);
end

end

