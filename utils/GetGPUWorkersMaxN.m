function n = GetGPUWorkersMaxN(volopt_)
%GETGPUWORKERSMAXN This function calculate the GPU workers max number
%depends on memory space using
% NOTE that each gpu processed volume has the same size: 256*256*25, which
% will allocate ~1GB graphics memory (uninclude worker self use)
% Input:
% Output:
%   - n: 1-by-1 positive integer, the maximum possible workers can be
%   allocated in this task

FOLD_CONST = 96;
SINGLE_BYTES = 4;

if ~exist("volopt_", "var")
    MEM_VOL_CONST = 2^30;
else
    MEM_VOL_CONST = volopt_.width * volopt_.height * volopt_.slices * ...
        SINGLE_BYTES * FOLD_CONST;
end

if ispc()
    MEM_WORKER_CONST = 800*1024*1024; % bytes
else
    MEM_WORKER_CONST = 1100*1024*1024; %bytes
end
MEM_SECURATY_RATIO = 0.85;
GPU_SECURATY_RATIO = 1;

n_gpu = gpuDeviceCount;

if n_gpu == 0
    % no available gpu
    n = 0;
    return;
end

gpu_memory = zeros(1, n_gpu);
for k = 1:n_gpu
    c = gpuDevice(k);
    gpu_memory(k) = c.TotalMemory;
end
minmem_tot = min(gpu_memory);


n = round(n_gpu*minmem_tot*MEM_SECURATY_RATIO/...
    (MEM_VOL_CONST+MEM_WORKER_CONST));

% gpu number can't be bigger than cpu number
n = round(min(n, GetCPUWorkersMaxN([],[]))*GPU_SECURATY_RATIO);

if n > n_gpu
    n = n - mod(n, n_gpu);
end

end