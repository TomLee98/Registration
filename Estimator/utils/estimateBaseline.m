function bl = estimateBaseline(F, opts)
%ESTIMATEBASELINE This function estimates the baseline of F by given
% estimation options
% Input:
%   - F: s-by-t double array, row for components, colume for time step
%   - opts: 1-by-1 sigopt object
% Output:
%   - bl: s-by-t double array, row for components, colume for time step
% 
% see also: movrank1

arguments
    F       (:,:)   double  {mustBePositive}
    opts    (1,1)   sigopt
end

switch opts.BaselineModel
    case "MovingQuantile"
        args = struct("q",      opts.Options.MinQuantileValue, ...
                      "auto",   opts.Options.AutoQuantile, ...
                      "win",    opts.Options.WindowSize);
        bl = estimateBaseline_MQ(F, args);
    case "MixedExponential"
        args = struct("q",      opts.Options.MinQuantileValue, ...
                      "auto",   opts.Options.AutoRegress, ...
                      "order",  opts.Options.RegressOrder);
        bl = estimateBaseline_ME(F, args);
    otherwise
end

end


function bl = estimateBaseline_MQ(F, args)

Estimator.auto_parpool("on");

bl = nan(size(F), "like", F);

if args.auto == true
    w = args.win;
    parfor k = 1:size(F,2)
        q = auto_quantile(F(k, :));
        bl(k, :) = movrank1(F(k, :), q, w);
    end
else
    if args.win >= size(F,2)
        bl = quantile(F, 2)*ones(1, size(F,2));
    else
        q = args.q;
        w = args.win;
        parfor k = 1:size(F,2)
            bl(k, :) = movrank1(F(k, :), q, w);
        end
    end
end

Estimator.auto_parpool("off");

end

function [qut, val] = auto_quantile(F)
%auto_prctile - Estimate Background Percentile.
% input:
%   - F: n-by-1 or 1-by-n double vector
% output:
%   - qut: the background quantile estimation
%   - val: the value compare with the quantile estimation
%   
% see also kde

[~, meshs, density, cdfs] = kde(F);
[~, argMax] = max(density);
qut = cdfs(argMax)*100;
val = meshs(argMax);

if isnan(qut)
    warning("NaN percentile computed. Reverting to median.")
    qut = 50;
    val = median(F);
end

end


function [bandWidth, meshs, density, cdfs] = kde(data, N)
%KDE - Kernal Density Estimate.
% This function for estimating smooth cdf by using ksdensity
% input:
%   - data: 1-D numeric vector
%   - N: bin counts
% output:
%   - bandWidth: the band width came from ksdensity
%   - meshs: the bin center vector
%   - density: the probability density
%   - cdfs: the cumulative probability density
%
%   see also ksdensity

% Parameters to set up the mesh on which to calculate
if ~exist("N","var")
    N = 2^12;
else
    N = max(N, 1024);
    N = fix(2^ceil(log2(N)));
end

[density,meshs,bandWidth] = ksdensity(data,[],"NumPoints",N);

cdfs = cumsum(density)*(meshs(2)-meshs(1));
end


function bl = estimateBaseline_ME(F, args)

end


