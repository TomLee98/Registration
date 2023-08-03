function [code,info] = GetHardwareInfo()
%GETHARDWAREINFO This function for detecting the hardware information
% Input:
% Output:
%   - code: a struct with 'cpu' and 'gpu', each item is a string
%   - info: a struct with 'cpu' and 'gpu', with the hardware information 
%     details

% Version 1.0.0
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
end

