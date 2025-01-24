% This script reads all volumes from ./raw_files and generates the random
% slices with the same resolution on XYZ
% Note: this script must be called in command window
%%
clear;
clc;

%% set parameters
resample_ratio = 0.5;           % random resample ratio at each dimension of volume
auto_resize = true;             % auto resize flag, true (strongly recommended) as default
suffix = "im";

%%
pathToThis = fileparts(mfilename('fullpath'));
src_folder = [pathToThis, '\raw_files'];
dst_folder = [pathToThis, '\trainset'];
fprintf("[source] %s\n[destination] %s\n", src_folder, dst_folder);

list = string({dir(src_folder).name});

istiff = list.endsWith(".tif","IgnoreCase",true);

if any(istiff)
    tiff_files = string(src_folder) + "\" + list(istiff);
    for k = 1:numel(tiff_files)
        %% load image data and information
        [V, vsiz, dV] = load(tiff_files(k));

        %% resize for cubic voxel
        if auto_resize, [V, vsiz, dV] = resize(V, vsiz, dV); end

        %% export random slices at each dimension
        rand_export(V, dst_folder, tiff_files(k), suffix, resample_ratio);
        if mod(k, 10) == 0, fprintf("%d volumes processed.\n", k); end
    end
end


%%

function [V, vsiz, dV] = load(file)
info = imfinfo(file); info = info(1);

vsiz = [info.Width, ...
        info.Height, ...
        str2double(string(info.ImageDescription).extractBetween("slices=",newline))];

dtype = sprintf("uint%d", info.BitDepth);

V = zeros(vsiz([2,1,3]), dtype);

dV = [1/info.XResolution, ...
    1/info.YResolution, ...
    str2double(string(info.ImageDescription).extractBetween("spacing=",newline))];

for s = 1:vsiz(3), V(:,:,s) = imread(file, s); end      % 2X faster than tiffReadVolume

end

function [V, vsiz, dV] = resize(V, vsiz, dV)
% padding to square
wp = max(vsiz(1), vsiz(2));
if vsiz(1) == vsiz(2)
    % already square image
    dw = [0, 0];
    dh = [0, 0];
else
    dw = [floor((wp - vsiz(2))/2), ceil((wp - vsiz(2))/2)];
    dh = [floor((wp - vsiz(1))/2), ceil((wp - vsiz(1))/2)];
end
V = padarray(V, [dh(1), dw(1)], "replicate","pre");
V = padarray(V, [dh(2), dw(2)], "replicate","post");

% padding some slices for cubic
% !!! assume that object length is larger on XY rather than Z
rV = uint8(rescale(V,0,255));
T = graythresh(rV);
mc_z = reshape(sum(rV-T*255,[1,2]),[],1);
mc_z = round(sum(mc_z.*(0.5:1:vsiz(3)-0.5)')./sum(mc_z));
z_bdsiz = round(wp/dV(3)*dV(1));
ds = z_bdsiz - vsiz(3);
ds = [floor(ds/2), ceil(ds/2)];
if any(ds < 0)
    % remove some slices
    ds = [floor((z_bdsiz-1)/2), ceil((z_bdsiz-1)/2)];
    V = V(:,:,max(mc_z-ds(1),1):min(mc_z+ds(2),vsiz(3)));
else
    V = padarray(V, [0, 0, ds(1)], "replicate","pre");
    V = padarray(V, [0, 0, ds(2)], "replicate","post");
end

% resize to cubic
V = imresize3(V, wp*ones(1,3), "Method","linear");

vsiz = size(V);
dV = dV(1)*ones(1,3);
end

function rand_export(V, dest, file, sfix, rr)
[~, name, ~] = fileparts(file);
prefile = fullfile(dest, name);

[l, w, h] = size(V);

sl = randperm(l, max(round(l*rr),1));
sw = randperm(w, max(round(w*rr),1));
sh = randperm(h, max(round(h*rr),1));

for y = sl
    imwrite(squeeze(V(y,:,:)), prefile+"_XZ_"+string(y)+"_"+sfix+".png", "png");
end
for x = sw
    imwrite(squeeze(V(:,x,:))', prefile+"_YZ_"+string(x)+"_"+sfix+".png", "png");
end
for z = sh
    imwrite(V(:,:,z), prefile+"_XY_"+string(z)+"_"+sfix+".png", "png");
end

end