function mem = GetAvailableMemory()
%GETAVAILABLEMEMORY This function get the available memory size with unit
%Bytes
% Output:
%   - mem: the available memory, unit with byte

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

% SELECT PRESENT PLATFORM
if ispc()
    [~, sysview] = memory();
    mem = sysview.PhysicalMemory.Available;
elseif isunix()
    [~,w] = unix('free | grep Mem');
    stats = str2double(regexp(w,'[0-9]*','match'));
    mem = (stats(3)+stats(end))*1e3;    % bytes
else
    % what is the fuck?
    error('Your operation system is so coooool.');
end

end

