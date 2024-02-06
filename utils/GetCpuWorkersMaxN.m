function n = GetCpuWorkersMaxN(volopt_, regopt_)
%GETCPUWORKERSMAXN This function calculate the CPU workers max number
%depends on memory space using
% Input
%   - volopt_: 1-by-12 table, with volumes information
%   - regopt_: 1-by-1 regopt object, registration options
% Output:
%   - n: 1-by-1 positive integer, the maximum possible workers can be
%   allocated in this task

UINT16_BYTES = 2;
% left ?% avaiable memory for system running AND O(1) memory alloc
MEM_SECURATY_RATIO = 0.8;

FOLD_RATIO = 5;

n_cpu = feature('numCores');

if  n_cpu <= 2
    % single or double core(s) cpu, no running supported
    n = 0;
    return;
elseif (n_cpu > 2) && (n_cpu <= 16)
    CPU_SECURATY_RATIO = 1;
else
    CPU_SECURATY_RATIO = 0.85;
end

if isempty(volopt_) && isempty(regopt_)
    mem_per_volume = 0;

    % SELECT PRESENT PLATFORM
    if ispc()
        mem_per_worker = 800*1024*1024; % bytes
        n = fix(GetAvailableMemory()*MEM_SECURATY_RATIO...
            /(mem_per_volume*FOLD_RATIO+mem_per_worker));
    else
        mem_per_worker = 1100*1024*1024; %bytes
        n = fix(GetAvailableMemory()*MEM_SECURATY_RATIO...
            /(mem_per_volume*FOLD_RATIO+mem_per_worker));
    end
    n = min(n, round(n_cpu*MEM_SECURATY_RATIO));
elseif istable(volopt_) && all(size(volopt_)==[1,12]) ...
        && isa(regopt_, "regopt")
    if volopt_.channels == 1

    else
        % THE MAXSIZE OF MEMORY CONTAINS MODEL CAN BE
        % CALCULATE BY LINEAR SIMILARITY
        mem_per_volume = volopt_.width*volopt_.height*volopt_.slices ...
            *UINT16_BYTES;

        % SELECT PRESENT PLATFORM
        if ispc()
            mem_per_worker = 800*1024*1024; % bytes
            n = fix(GetAvailableMemory()*MEM_SECURATY_RATIO...
                /(mem_per_volume*FOLD_RATIO+mem_per_worker));
        else
            mem_per_worker = 1100*1024*1024; %bytes
            n = fix(GetAvailableMemory()*MEM_SECURATY_RATIO...
                /(mem_per_volume*FOLD_RATIO+mem_per_worker));
        end

        if n < 1
            warning("GetCPUWorkersMaxN:noEnoughMemory", ...
                "No available enough memory.");
        end
    end

    n = min(n, round(n_cpu*CPU_SECURATY_RATIO));
else
    throw(MException("GetCPUWorkersMaxN:invalidInputArgs", ...
        "Invalid input arguments."));
end

end
