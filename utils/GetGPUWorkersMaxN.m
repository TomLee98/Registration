function n = GetGPUWorkersMaxN(volopts)
%GETGPUWORKERSMAXN This function calculate the GPU workers max number
%depends on space using
% input:
%   - volopts: 1-by-12 table, with volumes information
% output:
%   - n: 1-by-1 positive integer, the maximum possible workers can be
%   allocated in this task

arguments
    volopts (1,12) table
end

assert(volopts.channels <= 2, "GetGPUWorkersMaxN:invalidChannelsNumber", ...
    "Unsupported registration mode.");

SINGLE_BYTES = 4;
FOLD_RATIO = 160;   % the linear estimation may be wrong
SECURATY_RATIO = 0.85;

% THE MAXSIZE OF GRAPHICS CARD MEMORY CONTAINS MODEL CAN BE
% CALCULATE BY LINEAR SIMILARITY

gpu_memory = zeros(1,gpuDeviceCount);
for k = 1:gpuDeviceCount
    c = gpuDevice(k);
    gpu_memory(k) = c.AvailableMemory;
end
minimal_aval_memory = min(gpu_memory);

mem_per_volume = volopts.width*volopts.height*volopts.slices ...
    *SINGLE_BYTES;

if ispc()
    mem_per_worker = 800*1024*1024; % bytes
    n = round(gpuDeviceCount*minimal_aval_memory*SECURATY_RATIO/...
        (mem_per_volume*FOLD_RATIO+mem_per_worker));
elseif isunix()
    mem_per_worker = 1100*1024*1024; %bytes
    n = round(gpuDeviceCount*minimal_aval_memory*SECURATY_RATIO/...
        (mem_per_volume*FOLD_RATIO+mem_per_worker));
else
    n = 1;  % ? operation system
end

% for even gn to balance
if n < 1 || mem_per_volume > minimal_aval_memory
    throw(MException("GetGPUWorkersMaxN:noEnoughMemory", ...
        "No available enough memory."));
else
    if n > gpuDeviceCount
        % balance different GPU loading
        n = n - mod(n,gpuDeviceCount);
    end
end
end