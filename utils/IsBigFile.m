function res = IsBigFile(mem_used, frameratio, threshold)
%ISBIGFILE This function adjust whether loader allocate many memory, compare
%with current memory size and threshold (percentage)
% Input:
%   - mem_used: the raw image file used memory, bytes
%   - frameratio: the ratio of processing frames and total frames
%   - threshold: the percentage threshold of memory, 12.5 as default
% Output:
%   - res: bool, big file flag

% Version 1.0.0
%   *** support auto threshold
%
% Version 1.0.1
%   *** use memory used replace only disk file analysis

% Copyright (c) 2022-2023, Weihan Li

arguments
    mem_used (1,1) double {mustBeInteger, mustBePositive, mustBeFinite}
    frameratio (1,1) double {mustBeInRange(frameratio, 0, 1)} = 1;
    threshold (1,1) double {mustBeGreaterThanOrEqual(threshold, 0), ...
        mustBeLessThanOrEqual(threshold, 100)} = 12.5;
end

% read the available memory and file size
memory_available = GetAvailableMemory();    % Bytes

% auto threshold
if threshold==0 || threshold==100
    cpu_core_n = feature('numCores');
    threshold = (0.9/5.5 - (cpu_core_n*1024^3/memory_available)/5.5)*100;
end

if mem_used*frameratio/memory_available > threshold/100
    res = true;
else
    res =false;
end

end

