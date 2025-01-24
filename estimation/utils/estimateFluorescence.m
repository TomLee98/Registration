function F = estimateFluorescence(mov, mask, opts, comps, fc)
%ESTIMATEFLUORESCENCE This function estimate the coarse fluorescence
% Input:
%   - mov: 1-by-1 regmov object
%   - mask: m-by-n-by-p nonnegtive integer array, label array, 0 for
%           background, 1,2,...,etc indicate component label
%   - opts: 1-by-1 sigopt object, extraction options
%   - comps: c-by-1 positive integer vector, 
%   - fc: 1-by-1 functional channel indicator, could be "r","g" or "b"
% Output:
%   - F: c-by-t double matrix, #components = c, #timesteps = t

arguments
    mov     (1,1)   regmov
    mask    (:,:,:)         {mustBeNumeric}
    opts    (1,1)   sigopt
    comps   (1,:)   double  {mustBePositive, mustBeInteger}
    fc      (1,1)   string  {mustBeMember(fc, ["r","g","b"])} = "g"
end

F = nan(numel(comps), mov.MetaData.frames);

% start parpool
Estimator.auto_parpool("on");

fc = (fc == mov.MetaData.cOrder);
fn = mov.MetaData.frames;
bkg = opts.Options.Background;

keropt = struct("auto", opts.Options.AutoKernel, ...
                "kernel", opts.Options.Kernel);
kers = calc_kernels(mask, comps, keropt);

% some par-worker process ending with function handle release, but shared
% workers will lost it
% so we need to pre-load data to memory
mov_data = mov.Movie;

parfor n = 1:numel(comps)
    % for each worker, use pre-calculated kernel to estimate signal
    fpar = nan(1, fn);
    [kr, bdbox] = kers{n}{:};       % deal the arguments
    [rb, cb, hb] = bdbox2ds(bdbox);

    for t = 1:fn
        % model the fluorescence with constant expected camera background
        mov_cr = mov_data(rb, cb, fc, hb, t) - bkg; %#ok<PFBNS>
        mov_cr = cast(squeeze(mov_cr), "double");

        % weighted intensity, R^n->R
        fpar(t) = sum(mov_cr.*kr, "all");
    end

    F(n, :) = fpar;
end

Estimator.auto_parpool("off");
end

function kers = calc_kernels(mask, comps, keropt)
kers = cell(size(comps));

if keropt.auto == true
    % use volume coefficient of variation to estimate the best kernel
    stats = regionprops3(mask, "Volume");
    stats(stats.Volume==0, :) = [];
    rho = std(stats.Volume)/mean(stats.Volume);
    if rho <= 0.3
        kernel = "uniform";
    elseif rho > 0.3 && rho <= 0.6
        kernel = "gaussian";
    else
        kernel = "log";
    end
else
    kernel = keropt.kernel;
end

parfor n = 1:numel(comps)
    bw_mask = (mask==comps(n));
    stats = regionprops3(bw_mask, "BoundingBox");
    % generate kernel which shape as BoundingBox
    bdbox = stats.BoundingBox; % (y,x,z)
    if size(bdbox, 1)~=1
        throw(MException("estimateFluorescence:invalidROIs", ...
            "There are disconnected Domain ROIs with one label[%d].", n));
    end

    switch kernel
        case "uniform"
            kr = fspecial3("average", bdbox([5,4,6]));
        case {'gaussian', 'log'}
            kr = fspecial3(kernel, bdbox([5,4,6]), bdbox([5,4,6])/4);
        otherwise
            % average as default
            kr = fspecial3("average", bdbox([5,4,6]));
    end

    % masked kr
    [rb, cb, hb] = bdbox2ds(bdbox);
    kr_mask = bw_mask(rb, cb, hb);
    kr = kr.*kr_mask;

    % renormalization kernel, s.t. weights sum equals to 1
    kr = kr/sum(kr, "all");

    kers{n} = {kr, bdbox};
end

end

function [rb, cb, hb] = bdbox2ds(bd)
rb = round(bd(2):bd(2)+bd(5)-1);
cb = round(bd(1):bd(1)+bd(4)-1);
hb = round(bd(3):bd(3)+bd(6)-1);
end