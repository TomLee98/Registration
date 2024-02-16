function [n, src] = ParseWorkersNumber(regopt_)
% PARSEWORKERSNUMBER This function gets the maximum workers number with 
% memory size as constrain by parse the passed regopt object
% Input:
%   - regopt_: 1-by-1 regopt object
% Output:
%   - n: 1-by-1 nonnegtive integer, the workers number
%   - src: 1-by-1 string, the running resource, could be "cpu" or "cpu|gpu"
%
% see also: regopt, GetCPUWorkersMaxN, GetGPUWorkersMaxN

arguments
    regopt_ (1,1)   regopt
end

if regopt_.Algorithm == "MANREG"
    n = 1;
    src = "cpu";
else
    switch regopt_.Mode
        case "global"
            n = GetCPUWorkersMaxN([], []);
            src = "cpu";
        case "local"
            if regopt_.SubAlgorithm == "usual"
                % local & usual only support cpu current version(R2023b)
                n = GetCPUWorkersMaxN([], []);
                src = "cpu";
            else
                switch regopt_.Options.Hardware
                    case "cpu"
                        n = GetCPUWorkersMaxN([], []);
                    case "cpu|gpu"
                        n = GetGPUWorkersMaxN();
                    otherwise
                end
                src = regopt_.Options.Hardware;
            end
        otherwise
    end
end
end
