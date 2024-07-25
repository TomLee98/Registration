function tasks_ = parse_man(movsrc_, movtmpl_, regfrs_, volopt_, regopt_, batchsz_)
%PARSE_MAN This function parse the manual registration options and 
% generate tasks, where task can be passed to RegisterWorker object
% Input:
%   - movsrc_:  1-by-1 regmov object, placeholder inuse
%   - movtmpl_: 1-by-1 regtmpl object, the registration template
%   - regfrs_:  1-by-n positive integer array, the volumes absolute indices
%   - volopt_:  1-by-12 table, the volume properties, placeholder inuse
%   - regopt_:  1-by-1 regopt, the registration options
%   - batchsz_: 1-by-1 positive integer, the running batch size, placeholder inuse
% Output:
%   - tasks_:   1-by-k mQueue, the tasks queue (job) for running
%
% see also task, regtmpl, regopt

arguments
    movsrc_     (1,1)   regmov %#ok<INUSA>
    movtmpl_    (1,1)   regtmpl
    regfrs_     (1,1)   double {mustBePositive, mustBeInteger}
    volopt_     (1,12)  table %#ok<INUSA>
    regopt_     (1,1)   regopt
    batchsz_    (1,1)   double {mustBePositive, mustBeInteger} %#ok<INUSA>
end

tasks_ = mQueue();

% send the same template and registration options
tk = task(regfrs_, movtmpl_, regopt_);

tasks_.enqueue(tk);
end

