function parobj = OpenParpool(n)
%OPENPARPOOL This function open the parpool and return the parpool object, 
% note that this function uses cluster: Reg3D_Server as configuration
% Input
%   - n: positive integer, indicating the parpool size
% Output
%   - parobj: parpool object

% using default profile
pcl = parcluster(parallel.defaultProfile);
% make sure your cpu supports 'hyper-threads' technology
% pcl.NumThreads = 2;
% get present parpool
parobj = gcp("nocreate");

if isempty(parobj)
    parobj = parpool(pcl, [1,n], 'SpmdEnabled',false);
elseif parobj.NumWorkers ~= n
    % restart parpool
    delete(gcp);
    parobj = parpool(pcl, [1,n], 'SpmdEnabled',false);
end
end

