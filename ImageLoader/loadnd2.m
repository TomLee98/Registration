function [info, data] = loadnd2(file, tspan)
%LOADND2 This function load nd2 files with tspan
% Input:
%   - file: 1-by-1 string, the data full file path
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
% Output:
%   - info: struct with {opts, rt}, where opts is 1-by-12 metadata table,
%     rt is n-by-1 nonnegtive double time array
%
% if date is not needed but info, use
%       info = loadnd2(file)
% or you want the data, use
%       [info, data] = loadnd2(file, tspan)
% where you can pass the time span as tspan for continous time cut

arguments
    file    (1,1)   string
    tspan   (1,:)   double = []
end

nargoutchk(1, 2);

file = char(file);

% we need to use the bfmatlab package
status = bfCheckJavaPath(1);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);
bfInitLogging();
r = bfGetReader(file, 0);
info_bf = r.getMetadataStore();

% extract the information
width = info_bf.getPixelsSizeX(0).getValue();
height = info_bf.getPixelsSizeY(0).getValue();
channels = info_bf.getPixelsSizeC(0).getValue();
slices = info_bf.getPixelsSizeZ(0).getValue();
frames = info_bf.getPixelsSizeT(0).getValue();
images = channels*slices*frames;
xRes = double(info_bf.getPixelsPhysicalSizeX(0).value);
yRes = double(info_bf.getPixelsPhysicalSizeY(0).value);
if slices ~= 1
    zRes = double(info_bf.getPixelsPhysicalSizeZ(0).value);
else
    zRes = 0;   % only one plane, no resolution at z direction
end
dimOrder = string(split(info_bf.getPixelsDimensionOrder(0).getValue(),""));
dimOrder = dimOrder(strlength(dimOrder)>0)';
cOrder = getChannelOrder(info_bf);
dataType = string(info_bf.getPixelsType(0).getValue());
if dataType == "float"
    dataType = "single";
end

% generate opts
opts = table(width,...          % #pixel
            height,...          % #pixel
            channels,...        % #channels
            slices,...          % #slices
            frames,...          % #volumes
            images,...          % #images
            xRes,...            % μm/pixel
            yRes,...            % μm/pixel
            zRes,...            % μm/layer z
            dataType,...        % uint8, uint16, single
            dimOrder,...        % dimention order array
            cOrder);            % color channel order

rt = getRelativeTime(r, slices, channels, frames);

r.close();

info.opts = opts;
info.rt = rt;

if isempty(tspan), tspan = [1, opts.frames]; end

if nargout == 2
    if ispc()
        try
            % only PC support fast loading
            data = nd2Open3DVolume_reg(file, tspan);
        catch
            % no matter exception, call bfOpen3CVolume instead
            if opts.dimOrder(end) == "T"
                sspan = [(tspan(1)-1)*opts.slices*opts.channels+1, ...
                    tspan(2)*opts.slices*opts.channels];
                data = bfopen_reg(file, sspan);
            end
        end
    elseif isunix()
        data = bfOpen3DVolume_reg(file, tspan);
    end

    % reconstruct the image stack
    tmpopts = opts;
    tmpopts.frames = diff(tspan)+1;
    data = imreshape(data, tmpopts);
end

end

function mov = nd2Open3DVolume_reg(file, tspan)
% use nd2 library load total volumes
[FilePointer, ImagePointer, ImageReadOut] = ND2Open(file);
% mov is a c-by-1 cell, c for channels number
mov = ND2Read(FilePointer, ImagePointer, ImageReadOut, tspan(1):tspan(2));
ND2Close(FilePointer)
clear("FilePointer","ImagePointer","ImageReadOut");

for n = 2:numel(mov)
    mov{1} = cat(4, mov{1}, mov{n});
    mov{n} = [];    % free memory
end
mov = mov{1};
% permute as (Y,X,C,Z(T)): nd2 file fixed dimension order
mov = permute(mov, [1,2,4,3]);
end

function cOrder = getChannelOrder(info)
cc = info.getChannelCount(0);
C = zeros(cc,4,'uint8');
for k = 0:cc-1
    cinfo = info.getChannelColor(0,k);
    % some format there is no color channel information such as
    % tiff file exported from fiji
    if isempty(cinfo)
        cOrder = "";
        return;
    end
    C(k+1,:) = [cinfo.getRed(),cinfo.getGreen(),cinfo.getBlue(),...
        cinfo.getAlpha()];
end
C_na = C(:,1:3);  % cut the color map without alpha channel
% select the maximum value to represent the channel
[~,order] = max(C_na,[],2);
cOrder = ["r","g","b"];
cOrder = cOrder(order);
end

function rt = getRelativeTime(bf_reader, slices, channels, frames)
% get packaged meta data
store_meta_data = bf_reader.getMetadataStore();
rt = nan(frames, 1);
for k = 1:frames-1
    rt(k+1) = store_meta_data.getPlaneDeltaT(0,slices*channels*k-1).value();
end
rt = fillmissing(rt, "linear");
rt = round(rt, 3);
end

function mov = imreshape(mov, opts)
bindingArray = ["Y",    "X",    "C",    "Z",    "T";...
    "width","height","channels","slices","frames"];
[~,Loc] = ismember(opts.dimOrder,bindingArray(1,:));
mov = reshape(mov,...
    opts.(bindingArray(2,Loc(1))),...
    opts.(bindingArray(2,Loc(2))),...
    opts.(bindingArray(2,Loc(3))),...
    opts.(bindingArray(2,Loc(4))),...
    opts.(bindingArray(2,Loc(5))));
end