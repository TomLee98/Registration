function [u,v] = getBestDepose(s,w)
% if we have positive number s
% there we can find the positive number pair u,v
% s.t. min(w|s-uv|+(1-w)|u-v|),uv>=s,u>0,v>0,u,v is integer,0<=w<=1

arguments
    s (1,1) double {mustBePositive, mustBeInteger}
    w (1,1) double {mustBePositive, mustBeInRange(w,0,1)} = 0.5
end

z = sqrt(s);
u = [0;0;0];
v = [0;0;0];
% if ~exist('w','var') || ~isnumeric(w)||length(w)~=2 ...
%         || sum(w)~=1 || prod(w)<0
%     w = [0.5,0.5];
% end
w = [w, 1-w];
if z > 0
    for k = 1:3
        u(k) = fix(z)-2+k;
        if mod(s/u(k),1)~=0
            v(k) = fix(s/u(k))+1;
        else
            v(k) = s/u(k);
        end
    end
    % L determined by w(1): cover weight, and w(2): square weight
    L = w(1)*abs(u.*v-s)+w(2)*abs(u-v);
    p = find(L==min(L));
    p = p(randi(length(p),1));
end
u = u(p);
v = v(p);
end

