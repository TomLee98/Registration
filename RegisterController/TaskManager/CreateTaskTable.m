function T = CreateTaskTable(n)
%CREATETASKTABLE This function help to create a task table
% Input:
%   - n: 1-by-1 nonnegtive integer, rows of task table
% Output:
%   - T: n-by-6 table, with field {user, status, progress, time_used, n_cpu, n_gpu}
% see also: TaskManager

arguments
    n   (1,1)   double {mustBeNonnegative, mustBeInteger} = 0
end

T = table('Size', [n, 6], ...
    'VariableTypes', {'string','string','string','string','double','double'}, ...
    'VariableNames',{'user','status','progress','time_used','n_cpu','n_gpu'});
end
