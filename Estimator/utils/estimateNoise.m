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

dff_vid = dff(valid_idx, :);
noise_vid = noise(valid_idx, :);

parfor n = 1:numel(valid_idx)
    dff_n = dff_vid(n, :);

    % step1: estimate the components statistics properties
    switch nm
        case "normal"
            [p, ~] = eatimateNormalNoise(dff_n);

            w_sz = estimateNormalNoiseFilterWindow(p, dff_n);
            dff_n_sm = smoothdata(dff_n, "loess", round(w_sz));
            noise_n = dff_n - dff_n_sm;
        otherwise
            % do nothing
            dff_n_sm = dff_n;
            noise_n = zeros(size(dff_n));
    end

    dff_vid(n, :) = dff_n_sm;
    noise_vid(n, :) = noise_n;
end

dff(valid_idx, :) = dff_vid;
noise(valid_idx, :) = noise_vid;

Estimator.auto_parpool("off");

end

% ===================== utils functions =====================
% ======================== Normal Distribution =======================
function [x, fval] = eatimateNormalNoise(f)
v_max = max(f); v_min = min(f);

% assume that w1 = 0.9, w2 = 0.1
w10 = 0.9; w20 = 0.1;
g0 =  w20*(v_max-v_min)/mean(f);
d0 = max(-v_min, 0.01);
vsig_est = min(var(f)-0.01, w20*(v_max^2-v_min^2)/(2*g0) - w20^2/g0^2*(v_max-v_min)^2);
s0 = sqrt(var(f) - vsig_est);

% init parameters
x0 = [w10; w20; s0; g0; d0];
lb = [0;   0;   0;  0;  0.01];
ub = [1;   1;   1; inf; 10];
Aeq = [1,1,0,0,0];
beq = 1;

fminopts = optimoptions('fmincon','Algorithm','interior-point', ...
    'SubproblemAlgorithm','cg', 'EnableFeasibilityMode',true, ...
    'MaxIterations',round(0.5*numel(f)), ...
    'MaxFunctionEvaluations',round(5*numel(f)));

[x, fval] = fmincon(@(x)opfun_normal(x, f), x0, [],[], Aeq,beq, lb,ub, [], fminopts);
end

function f = opfun_normal(x, data)
% x: w1, w2, s, g, d
% p(u) = w1*1/(sqrt(2*pi)*s)*exp(-u^2/(2*s^2)) + w2*1/(g*u)

p1d = data(data <= x(5));
p2d = data(data > x(5));

p = prod(x(1)./(sqrt(2*pi)*x(3))*exp(-p1d.^2./(2*x(3)^2))) * ...
    prod(x(2)./(x(4)*p2d+eps));

f = -log(p);       % negtive Log-likelihood
end

function v = estimateNormalNoiseFilterWindow(p, f)
% apply lowess filter on f and minimize the noise properties
% use fminbnd for single variable: window_size optimization

nnp = [p(1), p(3)];

v = fminbnd(@(x)opfun_normal_winsize(x, nnp, f), 3, round(numel(f)/2));

end

function v = opfun_normal_winsize(w, p, f)
% v as Mahalanobis distance
% w: the window size need to be optimal
% p: noise releated coefficient, [weight, sigma]
% f: raw data

% local LSE for gaussian error removing
f_sm = smoothdata(f, "loess", round(w));

err = f - f_sm;

err = reshape(err, [], 1);

% use reformed Jensen–Shannon divergence
v = err'*err/var(err) + (var(err)+p(2)^2)/(std(err)*p(2)) - 2;
end
