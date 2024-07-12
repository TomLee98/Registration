function [code,info] = GetHardwareInfo()
%GETHARDWAREINFO This function for detecting the hardware information
% Input:
% Output:
%   - code: a struct with 'cpu' and 'gpu', each item is a string
%   - info: a struct with 'cpu' and 'gpu', with the hardware information 
%     details

% Version 1.1.0
% Copyright (c) 2022-2023, Weihan Li

% =============== process CPU information ==================
cpuinfo = getCpuInfo();
if contains(lower(cpuinfo.Name),'intel')
    code.cpu = "Intel";
elseif contains(lower(cpuinfo.Name),'amd')
    code.cpu = "AMD";
else
    code.cpu = "Unknown";
end
info.cpu = string(fieldnames(cpuinfo))+ ...
    ": " + string((struct2cell(cpuinfo)));

% =============== process GPU information ==================
gpuinfo = getGpuInfo();
if isempty(gpuinfo)
    code.gpu = "No GPU";
    info.gpu = "";
else
    if contains(lower(gpuinfo.Name),'intel')
        code.gpu = "Intel";
    elseif contains(lower(gpuinfo.Name),'amd')
        code.gpu = "AMD";
    elseif contains(lower(gpuinfo.Name),'nvidia')
        code.gpu = "Nvidia";
    else
        code.gpu = "Unknown";
    end
    info.gpu = string(fieldnames(gpuinfo))+ ...
        ": " + string((struct2cell(gpuinfo)));
end

% =============== process Memory information ==================
if ispc()
    [userview, sysview] = memory();
    code.mem = "general";
    info.mem = sprintf("Total Phys Mem: %d MB\n" + ...
                       "Available Phys Mem: %d MB\n" + ...
                       "MATLAB Used: %d MB", ...
                       round(sysview.PhysicalMemory.Total/1024/1024), ...
                       round(sysview.PhysicalMemory.Available/1024/1024), ...
                       round(userview.MemUsedMATLAB/1024/1024));
elseif isunix()
    code.mem = "general";
    [~,w] = unix('free -m | grep Mem');
    % total    used    free    shared    buff/cache    available
    stats = str2double(regexp(w,'[0-9]*','match'));
    [~, id] = unix('whoami');
    [~,w] = unix(['ps -aux | grep MATLAB | grep ', id]);
    w = string(w).splitlines();
    % user  pid  cpu_sys  cpu_user  ?  mem  group  pty  st  rt  path
    w = w(1).split(" ");   
    w(w=="")=[];
    info.mem = sprintf("Total Phys Mem: %d MB\n" + ...
                       "Available Phys Mem: %d MB\n" + ...
                       "MATLAB Used: %d MB", ...
                       stats(1), ...
                       stats(end), ...
                       str2double(w(6))/1024);
else
    % what is the fuck?
    throw(MException("GetAvailableMemory:unknownOS", ...
        "Unknown operation system."));
end

end

