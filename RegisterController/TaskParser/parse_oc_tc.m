function tasks_ = parse_oc_tc(movsrc_, movtmpl_, regfrs_, volopt_, regopt_, batchsz_)
%PARSE_OC_TC This function parse the two channels and one channel  
% registration options and generate tasks, where task can be passed 
% to RegisterWorker object
% Input:
%   - movsrc_:  1-by-1 regmov object, placeholder inuse
%   - movtmpl_: 1-by-1 regtmpl object, the registration template
%   - regfrs_:  1-by-n positive integer array, the volumes absolute indices
%   - volopt_:  1-by-12 table, the volume properties, placeholder inuse
%   - regopt_:  1-by-1 regopt, the registration options
%   - batchsz_: 1-by-1 positive integer, the running batch size
% Output:
%   - tasks_:   1-by-k mQueue, the tasks queue (job) for running
%
% see also task, regtmpl, regopt

arguments
    movsrc_     (1,1)   regmov %#ok<INUSA>     % no use in OC/TC mode
    movtmpl_    (1,1)   regtmpl
    regfrs_     (1,:)   double {mustBePositive, mustBeInteger}
    volopt_     (1,12)  table %#ok<INUSA>
    regopt_     (1,1)   regopt
    batchsz_    (1,1)   double {mustBePositive, mustBeInteger}
end

tasks_ = mQueue();

% split registration frames equally
nfrs = numel(regfrs_);
ntasks = ceil(nfrs/batchsz_);

for k = 1:ntasks
    regfrs = regfrs_((k-1)*batchsz_+1:min(k*batchsz_, nfrs));

    % send the same template and registration options
    tk = task(regfrs, movtmpl_, regopt_);

    tasks_.enqueue(tk);
end

end
