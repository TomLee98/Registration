function tasks_ = parse_man(movtmpl_, regfrs_, volopt_, regopt_, batchsz_)
%PARSE_MAN This function parse registration algorithm and options to task
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
    regfrs_     (1,1)   double {mustBePositive, mustBeInteger}
    volopt_     (1,12)  table %#ok<INUSA>
    regopt_     (1,1)   regopt
    batchsz_    (1,1)   double {mustBePositive, mustBeInteger}
end

tasks_ = mQueue();
tk = task(regfrs_, movtmpl_, regopt_);
tasks_.enqueue(tk);
end

