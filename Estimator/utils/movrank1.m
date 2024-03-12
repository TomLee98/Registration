function [M] = movrank1(A,r,varargin)
%MOVRANK - Moving percentile value
% Y = movrank1(X,r,K) for a vector X and positive integer scalar K
% computes a centered moving percentile r by sliding a window of length K along
% X. Each element of Y is the local percentile of the corresponding values of
% X inside the window, with Y the same size as X. When K is even, the
% window is centered about the current and previous elements of X. The
% sliding window is truncated at the endpoints where there are fewer than
% K elements from X to fill the window. The r is a number in range [0,100].
% 
% Y = movrank1(X,r,[NB NF]) for a vector X and nonnegative integers NB and
% NF computes a moving percentile along the length of X, returning the local
% percentile of the previous NB elements, the current element, and the next
% NF elements of X.
% 
% Y = movrank1(...,'Endpoints',ENDPT) controls how the percentile is
% calculated at the endpoints of X, where there are not enough elements
% to fill the window. ENDPT can be either a scalar numeric or logical
% value or one of the following:
% 
%     'shrink'    - (default) compute the percentile over the number of
%                   elements of X that are inside the window, effectively
%                   reducing the window size to fit X at the endpoints.
%     'discard'   - compute the percentile only when the window is filled
%                   with elements of X, discarding partial endpoint
%                   calculations and their corresponding elements in Y.
%                   This truncates the output; for a vector X and window
%                   length K, Y has length LENGTH(X)-K+1.
%                         
% When ENDPT is a scalar numeric or logical value, the missing elements
% of X inside the window are replaced with that value and Y remains the
% same size as X.
% 
% Y = movrank1(...,'Method',MTD) controls what is the method for
% calculating the local quantile in the sliding window. MTD can be 
% one of the following:
%
%   'exact'         - (default) to compute by sorting as explained below. 
%   'approximate'   - to use an approximation algorithm based on 
%                     t-digests.
%
%   Quantiles are specified using cumulative probabilities, from 0 to 1.
%   For an N element vector X, QUANTILE computes quantiles as follows:
%      1) The sorted values in X are taken as the (0.5/N), (1.5/N),
%         ..., ((N-0.5)/N) quantiles.
%      2) Linear interpolation is used to compute quantiles for
%         probabilities between (0.5/N) and ((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to quantiles
%         for probabilities outside that range.
% 
% Example: Compute a 5-point centered moving percentile.
%     t = 1:10;
%     x = [4 8 6 -1 -2 -3 -1 3 4 5];
%     yc = movrank1(x,50,5);
%     plot(t,x,t,yc);
% 
% Example: Compute a 5-point trailing moving percentile.
%     t = 1:10;
%     x = [4 8 6 -1 -2 -3 -1 3 4 5];
%     yt = movrank1(x,50,[4 0]);
%     plot(t,x,t,yt);
% 
% Example: Compute a 5-point trailing moving percentile, ignoring the first 4
% window shifts that do not contain 5 input elements.
%     x = [4 8 6 -1 -2 -3 -1 3 4 5];
%     yd = movrank1(x,50,[4 0],'Endpoints','discard');

% Version 1.0, 2021/12/28, by W.Li

% Set the default parameters
defaultK = 3;
defaultEndPoints = 'shrink';
defaultMethod = 'exact';

% Set the default parameters
p = inputParser;
validA = @(x)validateattributes(x,{'numeric'},{'real','vector'});
validR = @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',100});
validK = @(x)validateattributes(x,{'numeric'},{'vector','integer','>=',0});
validEndPoints = @(x)(isequal(x,'shrink') || isequal(x,'discard') ...
    || validateattributes(x,{'numeric'},{'real','scalar'}) ...
    || validateattributes(x,"logical","scalar"));
validMethod = @(x)(isequal(x,'exact') || isequal(x,'approximate'));
addRequired(p,'A',validA);
addRequired(p,'r',validR);
addOptional(p,'k',defaultK,validK);
addParameter(p,'endpoints',defaultEndPoints,validEndPoints);
addParameter(p,'method',defaultMethod,validMethod);
parse(p,A,r,varargin{:});

if numel(p.Results.k) == 1
    if mod(p.Results.k,2) == 1
        % the case for odd number
        kb = (p.Results.k - 1)/2;
        kf = kb;
    else
        % the case for even number
        kb = p.Results.k/2;
        kf = kb - 1;
    end
elseif numel(p.Results.k) == 2
    kb = p.Results.k(1);
    kf = p.Results.k(2);
else
    error("Invalid window size format.");
end

M = zeros(size(A),'like',A);

% sliding window
% May need parallel for speed up
lenA = numel(A);

switch p.Results.endpoints
    case 'shrink'
        for n = 1:lenA
            M(n) = quantile(A(max(1,n-kb):min(n+kf,lenA)),r/100,"Method",p.Results.method);
        end
    case 'discard'
        for n = kb+1:lenA-kf
            M(n-kb) = quantile(A(n-kb:n+kf),r/100,"Method",p.Results.method);
        end
        M(lenA-kf:end) = [];
    otherwise
        if isnumeric(p.Results.endpoints) && isscalar(p.Results.endpoints)
            if isrow(A) == true
                A = [repmat(p.Results.endpoints,1,kb),A,repmat(p.Results.endpoints,1,kf)];
            else
                A = [repmat(p.Results.endpoints,kb,1);A;repmat(p.Results.endpoints,kf,1)];
            end
            for n = kb+1:lenA+kb
                M(n-kb) = quantile(A(n-kb:n+kf),r/100,"Method",p.Results.method);
            end
        else
            error("Invalid end points format.");
        end
end

end
