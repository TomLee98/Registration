function fout = detrend_df_f(F,varargin)
%DETREND_DF_F - Compute Delta F over F.
%   This function returns the detrend Delta F over F df_f
%
%   fout = detrend_df_f(signal)
%   fout = detrend_df_f(___,Name,Value)
%
%   See also getSignal

% input:
%   - F: N by k matrix, with N frames and k components (neurons)
%   - varargin:
%       - WindowSize : number of frames for computing running 
%                      quantile, 30 (default)
%       - AutoQuantile: flag for determining quantile automatically, 
%                       true (default)
%       - MinQuantile: quantile used to estimate the baseline 
%                       (values in [0,100]) used only if 'flag_auto' is 
%                       False, i.e. ignored by default, 8 (default)
%       - UseMF: the flag for signal pre process whether apply mean
%                filter, value could be true,false,'auto', false as default
%       - OnlyF: the flag for only extract the fluoscrence intensity, true
%                or false, false as default
%       - Background: the estimated constant camera background, 
%                       usually about 100, 0 for no background extraction(default)
% output:
%   - df_f: the computed calcium activity to the derivative of F

% Version 1.0, 2021/12/30, by W.Li
%   *** basic extraction funciton
% Version 1.1, 2023/11/26, by W.Li
%   *** background option
%   *** only fluoscrence option

if nargin < 1
    error("Zero parameter is not supported.");
end

% Set the default parameters
minQuantile = 8;
windowSize = 30;
autoQuantile = true;
useMF = false;
onlyF = false;

% generate input checker
p = inputParser;
valid_signal = @(x)validateattributes(x,{'numeric'},{'2d','finite','nonnan','real'});
valid_min_quantile = @(x)validateattributes(x,{'numeric'},{'scalar','>',0});
valid_WindowSize = @(x)validateattributes(x,{'numeric'},{'scalar','>=',0});
valid_background = @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',65535});
valid_onlyF = @(x)validateattributes(x,{'logical'},{'scalar'});
addRequired(p,'signal',valid_signal);
addParameter(p,'MinQuantile',minQuantile,valid_min_quantile);
addParameter(p,'WindowSize',windowSize,valid_WindowSize);
addParameter(p,'AutoQuantile',autoQuantile,@(x)(islogical(x) && isscalar(x)));
addParameter(p,'UseMF',useMF,@(x)(islogical(x)||isstring(x)||ischar(x)));
addParameter(p,'Background',0,valid_background);
addParameter(p,'OnlyF',onlyF,valid_onlyF);
parse(p,F,varargin{:});

% data pre process
F = mean_filter(F,p.Results.UseMF,'DutyRatio',0.25);
bkg = p.Results.Background;
if p.Results.AutoQuantile == true
    if p.Results.WindowSize == 0 || p.Results.WindowSize > size(F,1)
        if p.Results.OnlyF
            fout = F-bkg;
        else
            % calculatet the detrend deltaF over F
            [dataPrct,~] = df_percentile(F,1);
            Fd = diag(prctile(F,dataPrct));
            fout = (F-Fd')./(Fd'-bkg);
        end
    else
        if p.Results.OnlyF
            fout = F-bkg;
        else
            windowSize = p.Results.WindowSize;
            [dataPrct,~] = df_percentile(F(1:windowSize,:),1);
            Fd = zeros(size(F),'like',F);
            for k = 1:size(F,2)
                Fd(:,k) = movrank1(F(:,k),dataPrct(k),windowSize);
            end
            fout = (F-Fd)./(Fd-bkg);
        end
    end
else
    if p.Results.WindowSize == 0 || p.Results.WindowSize > size(F,1)
        if p.Results.OnlyF
            fout = F-bkg;
        else
            Fd = prctile(F,p.Results.MinQuantile);
            fout = (F-Fd)./(Fd-bkg);
        end
    else
        if p.Results.OnlyF
            fout = F-bkg;
        else
            minQuantile = p.Results.MinQuantile;
            windowSize = p.Results.WindowSize;
            Fd = zeros(size(F),'like',F);
            for k = 1:size(F,2)
                Fd(:,k) = movrank1(F(:,k),minQuantile,windowSize);
            end
            fout = (F-Fd)./(Fd-bkg);
        end
    end
end

end

function signal = mean_filter(signal,varargin)
%MEAN_FILTER - Advanced mean filter 
%   This function apply mean filter on signal or not
%   
%   signal = mean_filter(signal)
%   signal = mean_filter(signal,mode)
%   signal = mean_filter(signal,mode,window)
%   signal = mean_filter(___,Name,Value)
%
%   See also movmean

% input:
%   - signal: N by k matrix, with N frames and k components
%   - varargin:
%       - (optional) WindowSize: the sliding window size, 5 (default)
%       - (optional) Mode: the filter running mode, true for apply filter
%                    on each sample, false means not, 'auto' for
%                    automatically apply filter on some samples if they
%                    satisfy the smooth condition number is large enough.
%                    Where the condition number is calculated by otsu-ALG.
%       - DutyRatio: enable only if Mode is 'auto', which is the controller
%                    of the automatic algorithm
% output:
%   - signal: processed N by k matrix, with N frames and k components
%
% Version 1.0, 2021/12/30, by W.Li

if nargin < 1
    error("Invalid input parameters.");
end

defaultWindowSize = 5;
defaultMode = 'auto';
defaultDutyRatio = 0.2;

p = inputParser;
addRequired(p,'signal',...
    @(x)validateattributes(x,{'numeric'},{'nonempty','2d','real'}));
addOptional(p,'Mode',defaultMode,@(x)(islogical(x)||ischar(x)||isstring(x)));
addOptional(p,'WindowSize',defaultWindowSize,...
    @(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
addParameter(p,'DutyRatio',defaultDutyRatio,...
    @(x)validateattributes(x,{'numeric'},{'nonempty','real','scalar','>',0,'<',1}));
parse(p,signal,varargin{:});

if isequal(p.Results.Mode,'auto')
    % use otsu algorithm to find the signal threshold
    signalZip = mapminmax(signal',0,255)';
    for k = 1:size(signal,2)
        counts = histcounts(signalZip(:,k),256);
        T = otsuthresh(counts);
        % if signal duty ratio is high enough, we can use mean filter
        if T > p.Results.DutyRatio
            signal(:,k) = movmean(signal(:,k),p.Results.WindowSize,'omitnan');
        end
    end
elseif p.Results.Mode == true
    signal = movmean(signal,p.Results.WindowSize,1,'omitnan');
end

end


function [dataPercentile,val] = df_percentile(data,dim)
%DF_PERCENTILE - Estimate Background Percentile.
%   Extracting the percentile of the data where the mode occurs and its value.
%   Used to determine the filtering level for DF/F extraction. Note that
%   computation can be innacurate for short traces.
%
%   dataPercentile = df_percentile(data)
%   dataPercentile = df_percentile(___,dim)
%   [___,val] = df_percentile(data)
%   [___,val] = df_percentile(___,dim)
%
%   See also kde

% input:
%   - data: N by k matrix, with N frames and k components
%   - (optional) dim: the direction when calculating percentile 
% output:
%   - dataPercentile: the background percentile estimation
%   - val: the value compare with the percentile estimation
%
% Version 1.0, 2021/12/30, by W.Li

if ~ismatrix(data)
    error("Invalid data format.");
end

if ~exist("dim","var") || isempty(dim)
    [~,meshs,density,cdfs] = kde(data);
    [~,argMax] = max(density);
    dataPercentile = cdfs(argMax)*100;
    val = meshs(argMax);
    
    if isnan(dataPercentile)
        warning("NaN percentile computed. Reverting to median.")
        dataPercentile = 50;
        val = median(data);
    end
else
    if dim == 1
        dataPercentile = zeros(1,size(data,2));
        val = zeros(1,size(data,2));
        for k = 1:size(data,2)
            [dataPercentile(k),val(k)] = df_percentile(data(:,k));
        end
    elseif dim == 2
        dataPercentile = zeros(size(data,1),1);
        val = zeros(size(data,1),1);
        for k = 1:size(data,1)
            [dataPercentile(k),val(k)] = df_percentile(data(k,:));
        end
    else
        error("Invalid input matrix dimension.");
    end
end
end


function [bandWidth,meshs,density,cdfs] = kde(data,N)
%KDE - Kernal Density Estimate.
% This function for estimate smooth cdf by using ksdensity
%
%   [bandWidth,meshs,density,cdfs] = kde(data)
%   [bandWidth,meshs,density,cdfs] = kde(___,N)
%
%   See also ksdensity

% input:
%   - data: 1-D numeric vector
%   - N: bin counts
% output:
%   - bandWidth: the band width came from ksdensity
%   - meshs: the bin center vector
%   - density: the probability density
%   - cdfs: the cumulative probability density

if nargin < 1
    error("Invalid input parameters.");
end
if ~isvector(data)
    warning("Data is not 1-D array, with explicit conversion: reshape.");
    % reshape as column vector
    data = reshape(data,[],1);
end

% Parameters to set up the mesh on which to calculate
if ~exist("N","var")
    N = 2^12;
else
    N = fix(2^ceil(log2(N)));
end

[density,meshs,bandWidth] = ksdensity(data,[],"NumPoints",N);

cdfs = cumsum(density)*(meshs(2)-meshs(1));
end
