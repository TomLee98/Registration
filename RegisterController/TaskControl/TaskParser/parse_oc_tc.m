function tasks_ = parse_oc_tc(movtmpl_, regfrs_, volopt_, regopt_, batchsz_)
%PARSE_OC_TC This function parse registration algorithm and options to task
%queue, where task can be passed to RegisterWorker object
% Input:
%   - movtmpl_:
%   - regfrs_:
%   - volopt_:
%   - regopt_:
%   - batchsz_:
% Output:
%   - tasks_:
%
% see also task, regtmpl, regopt

arguments
    movtmpl_    (1,1)   regtmpl
    regfrs_     (1,:)   double {mustBePositive, mustBeInteger}
    volopt_     (1,1)   struct
    regopt_     (1,1)   regopt
    batchsz_    (1,1)   double {mustBePositive, mustBeInteger}
end





end
