function info = getGpuInfo()
%This function get gpu info
%   - info:
%       - Name: the GPU name, like 'NVIDIA Geforce RTX 3090'
%       - Clock: the GPU clock rate, like 1695 MHz
%       - AvailableMem: the GPU Available Memory
%       - NumProcessors: the GPU counts
%       - DriverVersion: the GPU Driver version
%       - ToolkitVersion: the GPU Toolkit Version

if gpuDeviceCount == 0
    info = [];
    return;
end

Info = gpuDevice();
info.Name = Info.Name;
info.Clock = [num2str(Info.ClockRateKHz/1000,'%d'),' MHz'];
info.AvailableMem = [num2str(round(Info.AvailableMemory/1024/1024),'%d'),' MB'];
info.NumProcessor = gpuDeviceCount;
info.DriverVersion = num2str(Info.DriverVersion,'%.4f');
info.ToolkitVersion = num2str(Info.ToolkitVersion,"%d");
end

