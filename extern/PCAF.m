function fd = PCAF(f,c)
% PCAF This function apply the PCA filter on f, which mixture with signal and
% noise, here we assume that signal power ratio is large enough
% Input:
%   - f: n*1 double vector need to be denoise
arguments
    f (:, 1) double;
    c (1, 1) double {mustBePositive, mustBeInteger, mustBeInRange(c,1,5)} = 2;
end

n = numel(f);
% split the f series
m = 2*ceil(sqrt(n/10));
sp_A = 1:m:n;
%shift sp for center the denoised range
A_shift = round((n-sp_A(end)+1)/2);
sp_A = sp_A + A_shift;
% generate the sample matrix with k rows and m columns
A = zeros(length(sp_A)-1, m);
for sn = 1:length(sp_A)-1
    A(sn, :) = f(sp_A(sn):sp_A(sn+1)-1)';
end

% using PCA and select the first three principle components
[coeff, score, ~] = pca(A,"Centered",false,...
    "NumComponents",c,"Algorithm","als");
A_denoised = score*coeff';

% generate the thinning for border estimation
% maximum difference between head and tail is one 
mm = round(m/2);
sp_B = 1:mm:n;
if sp_B(end)<n
    sp_B = [sp_B, sp_B(end)+mm];
end
B = zeros(length(sp_B)-1, mm);
for sn = 1:length(sp_B)-1
    if sn==length(sp_B)-1 && sp_B(end) > n
        B(sn, :) = [f(sp_B(sn):end)', nan(1,mm-numel(f(sp_B(sn):end)))];
    else
        B(sn, :) = f(sp_B(sn):sp_B(sn+1)-1)';
    end
end
% using PCA and select the first three principle components
[coeff, score, ~] = pca(B,"Centered",false,...
    "NumComponents",c,"Algorithm","als");
B_denoised = score*coeff';

% combine the splited f
A = reshape(A_denoised',1,[]);
B = reshape(B_denoised',1,[]);

fd = [B(1:A_shift)';
      A';
      B(end-(n-numel(A)-A_shift)+1:end)'];
end
