function n = ParseProtectedWorkersNumber(alg_, src_, nwp_)
%PARSEPROTECTEDWORKERSNUMBER This function get the protected workers number
% due to algorithm, resource and passed protected workers number
% Input:
%   - alg_: 1-by-1 string, algorithm identity
%   - src_: 1-by-1 string, the running resource, could be "cpu" or "cpu|gpu"
%   - nwp_: 1-by-2 nonnegtive integer, max [cpu, gpu] workers number
% Output:
%   - n: 1-by-1 nonnegtive integer, the protected workers number
% see also: regopt, TaskManager

arguments
    alg_    (1,1)   string
    src_    (1,1)   string  {mustBeMember(src_, ["cpu", "cpu|gpu"])}
    nwp_    (1,2)   double  {mustBeNonnegative, mustBeInteger}
end

switch alg_
    case "MANREG"
        n = 0;
    otherwise
        if src_ == "cpu"
            n = nwp_(1);
        elseif src_ == "cpu|gpu"
            n = nwp_(2);
        end
end
end
