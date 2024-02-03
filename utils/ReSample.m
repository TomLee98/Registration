function [series_rs, scale] = ReSample(series, scale)
%DOWNSAMPLING The function downsampling volume series for registration 
% speed up, using 'linear' as downsampling method
% Input:
%   - series: ndarray, dimension from 2 to 5
%   - scale: 0-1 for volume downsampling, 1-inf for upsampling

if ~exist("scale","var")
    scale = 256./max(size(series, [1, 2]));
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

