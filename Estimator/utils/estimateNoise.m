function [dff, noise] = estimateNoise(dff, opts)
%ESTIMATENOISE This function use EM algorithm to split two mixed
%RV, which come from normal distribution and power distribution (p=-1)
% and return parameters
% Input:
%   - dff:
%   - opts:
% Output:
%   - dff:
%   - noise:

arguments
    dff     (:,:)   double
    opts    (1,1)   sigopt
end

noise = nan(size(dff), "like", dff);
nm = opts.Options.NoiseModel;

Estimator.auto_parpool("on");

% using EM algorithm
% assume that noise ~ X distribution and signal ~ power distribution
% maximize the Log-likelihood for arguments
valid_idx = find(~isnan(dff(:,1)));

for n = 1:numel(valid_idx)
    vid = valid_idx(n);
    dff_n = dff(vid, :); %#ok<PFBNS>

    % step1: estimate the components statistics properties
    switch nm
        case "normal"
            [p, ~] = eatimateNormalNoise(dff_n);
        case "exponential"
            [p, ~] = eatimateExponentialNoise(dff_n);
        case "gamma"
            [p, ~] = eatimateGammaNoise(dff_n);
        otherwise
    end

    % step2: optimize the filter parameters s.t. properties matching good
    switch nm
        case "normal"
            [w_sz, ~] = estimateNormalNoiseFilterWindow(p, dff_n);
            dff_n_sm = smoothdata(dff_n, "loess", w_sz);
            noise_n = dff_n - dff_n_sm;
        case "exponential"

        case "gamma"

        otherwise
    end

    dff(vid, :) = dff_n_sm;
    noise(vid, :) = noise_n;
end


Estimator.auto_parpool("off");

end

% ===================== utils functions =====================
% ======================== Normal Distribution =======================
function [x, fval] = eatimateNormalNoise(f)
v_max = max(f); v_min = 0;

% assume that w1 = 0.9, w2 = 0.1
w10 = 0.9; w20 = 0.1;
g0 =  w20*(v_max-vmin)/mean(f);
s0 = sqrt(var(f) - w20*(v_max^2-v_min^2)/(2*g0));

% init parameters
x0 = [w10, w20, s0, g0];
lb = [0,   0,  0,   0];
ub = [1,   1,  inf, inf];
Aeq = [1,1,0,0];
beq = 1;

fminopts = optimoptions('fmincon','Algorithm','interior-point', ...
    'SubproblemAlgorithm','cg', 'EnableFeasibilityMode',true, ...
    'MaxIterations',round(0.5*numel(f)), ...
    'MaxFunctionEvaluations',round(5*numel(f)));

[x, fval] = fmincon(@(x)opfun_normal(x, f), x0, [],[], Aeq,beq, lb,ub, [], fminopts);
end

function f = opfun_normal(x, data)
% x: w1, w2, s, g
% p(u) = w1*1/(sqrt(2*pi)*s)*exp(-u^2/(2*s^2)) + w2*1/(g*u)

p = x(1)/(sqrt(2*pi)*x(3))*exp(-data.^2/(2*x(3)^2)) + ...
    x(2)*1/(x(4)*data);

f = -log(prod(p));       % negtive Log-likelihood
end

function v = estimateNormalNoiseFilterWindow(p, f)
% apply lowess filter on f and minimize the noise properties
% use fminbnd for single variable: window_size optimization

nnp = [p(1), p(3)];

v = fminbnd(@(x)opfun_normal_winsize(x, nnp, f), 3, numel(f)/2);

end

function v = opfun_normal_winsize(w, p, f)
% v as Mahalanobis distance
% w: the window size need to be optimal
% p: noise releated coefficient, [weight, sigma]
% f: raw data

% local LSE for gaussian error removing
f_sm = smoothdata(f, "loess", w);

err = f - f_sm;

err = reshape(err, 1, []);

v = sqrt(err'*err)/(p(1)*p(2));
end

% ================== Exponential Distribution ==================

function [x, fval] = eatimateExponentialNoise(f)

end

function f = opfun_exp(x, data)

end

% ================== Gamma Distribution ==================

function [x, fval] = eatimateGammaNoise(f)

end

function f = opfun_gamma(x, data)

end