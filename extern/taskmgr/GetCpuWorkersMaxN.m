function n = GetCpuWorkersMaxN(volopts, big_flag, longterm_flag)
%GETCPUWORKERSMAXN This function calculate the CPU workers max number
%depends on space using
% input:
%   - volopts: 1-by-12 table, with volumes information
%   - big_flag: 1-by-1 logical, indicate if open big-file mode
% output:
%   - n: 1-by-1 positive integer, the maximum possible workers can be
%   allocated in this task

arguments
    volopts (1,12) table
    big_flag (1,1) logical = false
    longterm_flag (1,1) logical = false
end

assert(volopts.channels <= 2, "GetCPUWorkersMaxN:invalidChannelsNumber", ...
    "Unsupported registration mode.");

if volopts.channels == 1

else
    UINT16_BYTES = 2;
    % left ?% avaiable memory for system running AND O(1) memory alloc
    SECURATY_RATIO = 0.8;

    if longterm_flag == false
        if big_flag == true
            FOLD_RATIO = 2+3;
        else
            % 2 for two channel(aligned,signal),2 for argument and temporary
            % variable (pyramid levels algorithm), 3 for tmpdata (coreg3D)
            FOLD_RATIO = 2*2+3;
        end
    else
        FOLD_RATIO = 5;
    end

    % THE MAXSIZE OF MEMORY CONTAINS MODEL CAN BE
    % CALCULATE BY LINEAR SIMILARITY
    mem_per_volume = volopts.width*volopts.height*volopts.slices ...
        *UINT16_BYTES;

    % SELECT PRESENT PLATFORM
    if ispc()
        mem_per_worker = 800*1024*1024; % bytes
        n = fix(GetAvailableMemory()*SECURATY_RATIO...
            /(mem_per_volume*FOLD_RATIO+mem_per_worker));
    elseif isunix()
        mem_per_worker = 1100*1024*1024; %bytes
        n = fix(GetAvailableMemory()*SECURATY_RATIO...
            /(mem_per_volume*FOLD_RATIO+mem_per_worker));
    else
        n = 1;  % ? operation system
    end

    if n < 1
        throw(MException("GetCPUWorkersMaxN:noEnoughMemory", ...
            "No available enough memory."));
    end

    [~, wkn] = GetRunningMode();
    n = min(n, wkn);
end

end
