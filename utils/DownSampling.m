function [series_ds, scale] = DownSampling(series, scale)
%DOWNSAMPLING The function downsampling volume series for registration 
% speed up, using 'nearest' as downsampling method
% Input:
%   - series: ndarray, dimension from 2 to 5
%   - scale: 0-1 for volume downsampling, 1-inf for upsampling
if ~exist("scale","var")
    scale = 256./max(size(series, [1, 2]));
end
ndim = ndims(series);

switch ndim
    case 2
        series_ds = imresize(series,scale,"nearest");
    case 3
        series_ds = imresize3(series,"Method","nearest",...
            "Scale",[scale,scale,1]);
    case 4
        series_tmp = imresize3(series(:,:,:,1),...
                "Method","nearest","Scale",[scale,scale,1]);
        series_ds = zeros([size(series_tmp, [1,2]),size(series,[3,4])],...
            "like",series);
        series_ds(:,:,:,1) = series_tmp;
        for n = 2:size(series, 4)
            series_ds(:,:,:,n) = imresize3(series(:,:,:,n),...
                "Method","nearest","Scale",[scale,scale,1]);
        end
    case 5
        series_tmp = imresize3(squeeze(series(:,:,1,:,1)),...
                "Method","nearest","Scale",[scale,scale,1]);
        series_ds = zeros([size(series_tmp, [1,2]), 2, size(series,[4,5])],...
            "like",series);
        for c = 1:2
            for n = 1:size(series,5)
                series_ds(:,:,c,:,n) = imresize3(squeeze(series(:,:,c,:,n)),...
                    "Method","nearest","Scale",[scale,scale,1]);
            end
        end
    otherwise
        warning("Unsupported downsampling.");
        series_ds = series;
end

end

