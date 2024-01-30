function [opts,rt,mptr,filename] = loadfile(filename, omitif, order_out, memmap, num_history)
%LOADFILE: This function for loading data (.ims,.tif,.nd2)
% input:
%   - filename:char array or string, the movie file name (with full path), can be empty
%   - omitif: the flag for omit tiff file information lost, false as default
%   - order_out: the output stacks order, ["X","Y","C","Z","T"] as default
%   - memmap: the flag for memory mapping enabled
%   - num_history: the number of movies history kept
% output:
%   - opts: some useful movie information, 1 X 12 table, with variables:
%           width   height  channels slices  frames images  xRes    yRes
%           zRes    bitDepth    dimOrder    cOrder
%   - mptr: the mpimg object, which is dynamic image file object
%   - rt: the relative camera time, n*1 array, seconds as unit
%   - filename: the file name, char array
%
%   see also: imload

% Copyright (c) 2022-2024, Weihan Li
% LOADFILE: 
% Verison: 1.0.0
%   *** basic loading function

arguments
    filename        string
    omitif      (1,1) logical = false
    order_out   (1,5) string = ["Y","X","C","Z","T"]
    memmap      (1,1) logical = true
    num_history (1,1) double {mustBeNonnegative, mustBeInteger} = 0
end

nargoutchk(1, 4);

opts = [];
mptr = [];
rt = [];

% check whether file is existing or not
if isempty(filename)
    % uigetfile seletion
    [file,path] = uigetfile(...
        {'*.ims;*.tif;*.nd2','volume files (*.ims,*.tif,*.nd2)'},...
        'volume data selector');
    if isequal(file,0) || isequal(path,0)
        return;
    else
        filename = fullfile(path,file);
    end
end

if nargout < 3
    [status, info] = imload(filename, true, omitif);
else
    [status, info, mov] = imload(filename, true, omitif);
end

if status == -1
    error("Loading failed.");
else
    opts = info.opts;
    rt = info.rt;

    if nargout >= 3
        % using order_out for transformed order
        mov = reorder(mov, opts.dimOrder, order_out);

        if memmap == true
            wbar = loadbar("memmory mapping");

            % tmpfolder = "E:\";  % debug
            % generate a new mapping file
            tmpfolder = mpimg.findtmpfolder(opts, num_history);
            mptr = mpimg(tmpfolder, [], mov, order_out);

            delete(wbar);
        else
            mptr = mov;
        end
            
        disp("Loading success.");
    end
end
end

function mov = reorder(mov, order_old, order_new)
[~, order] = ismember(order_new, order_old);
mov = permute(mov, order);
end