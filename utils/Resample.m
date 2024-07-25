function [scale, series_rs] = Resample(series, scale)
%DOWNSAMPLING The function downsampling volume series for registration 
% speed up, using 'linear' as downsampling method, series must with
% dimension order: X,Y(,C,Z,T)
% Input:
%   - series: ndarray, dimension from 2 to 5
%   - scale: 0-1 for volume downsampling, > 1for upsampling, inf for auto
%            sampling
% Output:
%   - scale: resampling scale
%   - series_rs: resampling movie series
%
% example 1:
% If you only need the auto resampling scale, use:
%   scale = ReSample(series) or scale = ReSample(series, inf)
%
% example 2:
% If you need the resampling series, use:
%   [scale, series_rs] = ReSample(series, scale)
%
% see also: imresize, imresize3

arguments
    series  (:,:,:,:,:)   uint16
    scale   (1,1)   double {mustBeNonnegative} = inf
end

EC50 = 256;
if isinf(scale)
    sz_max = max(size(series, [1, 2]));
    scale = EC50/(EC50 + sz_max);
end

if nargout == 1
    series_rs = [];
    return;
end

ndim = ndims(series);

switch ndim
    case 2
        series_rs = imresize(series,scale,"linear");
    case 3
        series_rs = imresize3(series,"Method","linear",...
            "Scale",[scale,scale,1]);
    case 4
        series_tmp = imresize3(series(:,:,:,1),...
                "Method","linear","Scale",[scale,scale,1]);
        series_rs = zeros([size(series_tmp, [1,2]),size(series,[3,4])],...
            "like",series);
        series_rs(:,:,:,1) = series_tmp;
        for n = 2:size(series, 4)
            series_rs(:,:,:,n) = imresize3(series(:,:,:,n),...
                "Method","linear","Scale",[scale,scale,1]);
        end
    case 5
        series_tmp = imresize3(squeeze(series(:,:,1,:,1)),...
                "Method","linear","Scale",[scale,scale,1]);
        series_rs = zeros([size(series_tmp, [1,2]), 2, size(series,[4,5])],...
            "like",series);
        for c = 1:2
            for n = 1:size(series,5)
                series_rs(:,:,c,:,n) = imresize3(squeeze(series(:,:,c,:,n)),...
                    "Method","linear","Scale",[scale,scale,1]);
            end
        end
    otherwise
        warning("Unsupported downsampling.");
        series_rs = series;
end

end

