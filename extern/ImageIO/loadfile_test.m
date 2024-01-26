function [opts,rt,memptr,filename] = loadfile_test(filename, tmpfolder, omitif, order_out)
%LOADFILE: This function for loading data (.ims,.tif,.nd2)
% input:
%   - filename:char array or string, the movie file name (with full path), can be empty
%   - omitif: the flag for omit tiff file information lost, false as default
%   - order_out: the output stacks order, ["X","Y","C","Z","T"] as default
% output:
%   - opts: some useful movie information, 1 X 12 table, with variables:
%           width   height  channels slices  frames images  xRes    yRes
%           zRes    bitDepth    dimOrder    cOrder
%   - memptr: the memmapper object, which is file manage interface
%   - rt: the relative camera time, n*1 array, seconds as unit
%   - filename: the file name, char array
%
%   see also: imload

% Copyright (c) 2022-2024, Weihan Li
% LOADFILE: 
% Verison: 1.0.0
%   *** basic loading function

arguments
    filename string
    tmpfolder (1,1) string {mustBeFolder} = "E:\si lab\Matlab Projects\Registration";
    omitif (1,1) logical = false
    order_out (1,5) string = ["X","Y","C","Z","T"]
end

nargoutchk(1, 4);

opts = [];
memptr = [];
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

        memptr = mpimg(tmpfolder, mov, order_out, true);

        disp("Loading success.");
    end
end


    function mov = reorder(mov, order_old, order_new)
        [~, order] = ismember(order_new, order_old);
        mov = permute(mov, order);
    end
end