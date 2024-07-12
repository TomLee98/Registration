function tasks_ = parse_lt(movsrc_, movtmpl_, regfrs_, volopt_, regopt_, batchsz_)
%PARSE_LT This function parse the long term registration options and
%generate tasks, where task can be passed to RegisterWorker object
% Input:
%   - movsrc_:  1-by-1 regmov object, registration chain needed
%   - movtmpl_: 1-by-1 regtmpl object, the registration template
%   - regfrs_:  1-by-n positive integer array, the volumes absolute indices
%   - volopt_:  1-by-12 table, the volume properties, placeholder inuse
%   - regopt_:  1-by-1 regopt, the registration options
%   - batchsz_: 1-by-1 positive integer, the running batch size
% Output:
%   - tasks_:   1-by-k mQueue, the tasks queue (job) for running
%
% see also task, regtmpl, regopt

% Note that the registration options in tasks are calculated in this scope
% for step1(registration chain generation). So that ltreg only do 
% step 2 (each volume fine tuning)

arguments
    movsrc_     (1,1)   regmov  
    movtmpl_    (1,1)   regtmpl
    regfrs_     (1,:)   double {mustBePositive, mustBeInteger}
    volopt_     (1,12)  table %#ok<INUSA>
    regopt_     (1,1)   regopt
    batchsz_    (1,1)   double {mustBePositive, mustBeInteger}
end

% generate a registration chain (common using, handle-likely)
rc = regchain(movsrc_, movtmpl_, regopt_.Options);

rc.generate();

% update movtmpl (possible)
movtmpl_ = rc.Template;

% modify the registration options
regopt_ = set(regopt_, "RegChain", rc);
regopt_ = set(regopt_, "Keyframes", rc.ChainInfo.data.kfidx);

% allocate the tasks
tasks_ = mQueue();

% split registration frames equally
nfrs = numel(regfrs_);
ntasks = ceil(nfrs/batchsz_);

for k = 1:ntasks
    regfrs = regfrs_((k-1)*batchsz_+1:min(k*batchsz_, nfrs));

    % send the same template and modofied registration
    tk = task(regfrs, movtmpl_, regopt_);

    tasks_.enqueue(tk);
end

end

