function CloseParpool(parobj)
%CLOSEPARPOOL This function closes a parpool and removes failed job records, 
% note that this function uses cluster: Reg3D_Server as configuration
% Input:
%   - parobj: parpool object
% Output: None

% close parpool
delete(parobj);

% remove the possible latest failed jobs
pcl = parcluster(parallel.defaultProfile);
delete(pcl.Jobs);
end