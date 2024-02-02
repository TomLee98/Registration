function [rm, nc, ng] = GetRunningMode(disp_flag)
%GETRUNNINGMODE this function detects the hardware and output the best section

arguments
    disp_flag (1,1) logical = false
end

cpu_cores_n=feature('numCores');
if cpu_cores_n <= 8
    nc = cpu_cores_n;
else
    % note that about 85% of total workers is experimental value
    nc = round(0.85*cpu_cores_n);
end
ng = gpuDeviceCount('available');

if ng > 0
    if ng == 1
        rm = 'gpu';
    else
        rm = 'multi-gpu';
    end
else
    if cpu_cores_n > 1
        rm = 'multi-cpu';
    else
        % single core for debugging
        rm = 'cpu';
    end
end

if disp_flag == true
    fprintf('-> available cpu core: %d \n-> available gpu: %d\n',...
        cpu_cores_n, ng);
    fprintf('-> maximum stable workers number: %d\n',nc);
end

% Comment line for debugging
% rm = 'cpu';
end

