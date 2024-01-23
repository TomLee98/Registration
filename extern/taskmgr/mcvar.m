function [mp, nld, t0] = mcvar(vols, ds, bkg, w)
%MCVAR This function calculates the mass center variation in movie
% For an information on 'translation' estimation source using
% input:
%   - vols: m-by-n-by-p-by-t 4-dimension array, contains raw data
%   - ds: the xy plane downsampling ratio, for fast mass center estimated,
%         2 as default 
%   - bkg: 1-by-1 double,  to determine which pixel will be counted as
%          foreground
%   - w: 1-by-1 logical, if true, the intensity of pixel will be
%          considered
% output:
%   - mp: 1-by-3 array, the mass center mean on x,y,z
%   - nld: 1-by-1 nonnegtive, the nonlinear translation cost
%   - t0: 1-by-1 positive integer, indicate the best location on time line
%         which can be template for 'translation' transformation

arguments
    vols (:,:,:,:) uint16 {mustBeInteger, mustBeNonnegative}
    ds   (1,1) double {mustBeInteger, mustBePositive} = 2 
    bkg  (1,1) double {mustBeNonnegative} = 100
    w    (1,1) logical = false
end

% downsampling on XY plane
vols = vols(1:ds:end,1:ds:end,:,:);

ps = size(vols);
yps = ps(1); 
xps = ps(2); 
zps = ps(3); 
tps = ps(4);

MC = zeros(tps, 3);

if w == false
    vols = single(vols - 1.1*bkg);      % 1.1 for noise level coefficient
    vols(vols>0) = 1; 
else
    vols = single(vols - bkg);
end

[X,Y,Z] = meshgrid(1:xps, 1:yps, 1:zps);

for t = 1:tps
    Vol_Proj = vols(:,:,:,t).*X;
    MC(t,1) = mean(Vol_Proj(Vol_Proj>0), "all");
    Vol_Proj = vols(:,:,:,t).*Y; 
    MC(t,2) = mean(Vol_Proj(Vol_Proj>0), "all");
    Vol_Proj = vols(:,:,:,t).*Z; 
    MC(t,3) = mean(Vol_Proj(Vol_Proj>0), "all");
end

mp0 = mean(MC);
mp = [mean(MC(:,1:2))*ds, mean(MC(:,3))];

[t0, nld]= calc_best_loc_t(MC, mp0);

end

function [t0, nld] = calc_best_loc_t(MC, mp, dt)
% This function find the best location on time t0, which minimal the local
% lost function: f = dist(v(t0), v_mean) + std([v(t0-dt), ..., v(t0+dt)])
% where v is mass center of volume, v(t0): at t0

arguments
    MC (:,3) single
    mp (1,3) single
    dt (1,1) double {mustBeInteger, mustBeNonnegative} = 2
end

D = sqrt(sum((MC-mp).^2, 2));
V = nan(size(D));

if dt > 0
    for t = 1:numel(V)
        Dw = D(max(t-dt,1):min(t+dt,numel(V))); % windowed D(istance) vector
        V(t) = std(Dw);
    end
    [~, t0] = min(D + V);   % near mass center and stable enough
else
    [~, t0] = min(D);
end

% calculate nonlinear sp
NLD = D;
NLD(D<1) = 1;   % 1 pixel as limit
nld = sum(log(NLD));

end