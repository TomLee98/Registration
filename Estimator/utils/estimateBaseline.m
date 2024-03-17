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
    F       (:,:)   double
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
                      "grid",   opts.Options.GridStep, ...
                      "auto",   opts.Options.AutoRegress, ...
                      "order",  opts.Options.RegressOrder);
        bl = estimateBaseline_ME(F, args);
    otherwise
end

end


function bl = estimateBaseline_MQ(F, args)
% this function use move quantile algorithm to calculate baseline

Estimator.auto_parpool("on");

bl = nan(size(F), "like", F);

if args.auto == true
    w = args.win;
    parfor k = 1:size(F, 1)
        if all(~isnan(F(k, :)))
            q = auto_quantile(F(k, :));
            bl(k, :) = movrank1(F(k, :), q, w);
        end
    end
else
    if args.win >= size(F,2)
        bl = quantile(F, args.q/100, 2)*ones(1, size(F,2));
    else
        q = args.q;
        w = args.win;
        parfor k = 1:size(F, 1)
            if all(~isnan(F(k, :)))
                bl(k, :) = movrank1(F(k, :), q, w);
            end
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
% this function models the baseline as mixed exponential function as 
% fluorescence decay based on fluorescence protein life time

% extract the data lower than quantile threshold

Estimator.auto_parpool("on");

bl = nan(size(F), "like", F);

if args.auto == true

else
    q = args.q;
    g = args.grid;
    d = args.order;
    parfor k = 1:size(F, 1)
        tmpf = F(k, :);

        if all(~isnan(tmpf))
            % calculate piecewise fluorescence baseline
            N = ceil(numel(tmpf)/g);
            for gs = 1:N
                gvidx = (gs-1)*g+1:min(gs*g, numel(tmpf));
                tmpf_gs = tmpf(gvidx);
                vq = quantile(tmpf_gs, q/100);
                tmpf_gs(tmpf_gs > vq) = nan;
                tmpf(gvidx) = tmpf_gs;
            end

            % use makima interpolation to replace the nan value
            % tmpf = fillmissing(tmpf, "makima");

            tmpx = (0:numel(tmpf)-1);
            tmpxx = tmpx;

            tmpxx(isnan(tmpf)) = [];
            tmpf(isnan(tmpf)) = [];

            % fit the mixed exponential function by LSE
            s_fit = gen_mixexp(d);
            st0 = [repmat(mean(tmpf)/d, 1, d); zeros(1, d)];
            st0 = reshape(st0, 1, []);

            ft = fittype(s_fit);
            ftopts = fitoptions(ft);
            ftopts = fitoptions(ftopts, "Robust", "Bisquare", ...
                "Lower", zeros(1, 2*d), "StartPoint", st0);

            fitobj = fit(tmpxx', tmpf', ft, ftopts);

            bl(k, :) = feval(fitobj, tmpx);
        end
    end
end

Estimator.auto_parpool("off");
end

function s = gen_mixexp(order)
arguments
    order   (1,1)   {mustBeMember(order, 1:5)}  % avoid too high order overfitting
end

s = "a1*exp(-b1*x)";

for k = 2:order
    s = s + sprintf("+a%d*exp(-b%d*x)", k, k);
end

end
