function [status, path] = savefile(I,filepath,varargin)
%SAVETIFF This function save 4-D or 5-D matrix as selected file format,
% support tiff, mat(ndmatrix)
% input:
%   - I: the 4-D matrix, default order is X,Y,Z,T, single/uint16
%   - filepath: the path you want to put file in
%   - varargin:
%       - (optional) dimorder: the order of dimension, which is the arrangement of
%                              'X','Y','Z','T','C', and 'XYCZT' as default
%       - (optional) compression: specifies the compression to use when 
%                                 writing the OME-TIFF file, no compression
%                                 as default
%       - (optional) metadata: struct, allows to use a custom description
%                    in minimalMetadata, which append to the key:
%                    ImageDescription.
%       - (optional) filename: the default file name, string
% output:
%   - status: the return status, SAVE_SUCCESS = 0; SAVE_FAILED = -1;
%   - path: the file saving folder
%
%   see also: bfsave

% Copyright (c) 2021, Weihan Li
% SAVEFILE: 
% Version: 1.0.0
%   *** basic saving function
%   *** status indicate the function status
%   *** many export data format support

SAVE_SUCCESS = 0;
SAVE_FAILED = -1;

if nargin == 0
    error('A matrix is required');
end
assert(ndims(I)==4,"savefile:invalidDimension","I is not 4D array.");

dimOrder = 'XYCZT';
compr = '';
fname = "signal_aligned";

n = 1;
while n < numel(varargin)
    switch varargin{n}
        case 'dimorder'
            dimOrder = varargin{n+1};
        case 'compression'
            compr = varargin{n+1};
        case 'metadata'
            metadata = varargin{n+1};
        case 'filename'
            fname = varargin{n+1};
        otherwise
            warning('invalid parameters');
    end
    n = n + 2;
end

file_suppoted = {'*.tif','Tag Image File Format (*.tif)';...
                             '*.mat','MAT file (*.mat)';...
                             '*.avi','Audio Video Interleaved file (*.avi)';...
                             '*.mp4','MPEG-4 Part 14 (*.mp4)'};

if ~exist("filepath","var") || isempty(filepath) || ~exist(filepath,"dir")
    [file,path] = uiputfile(file_suppoted,'Save',fname);
else
    if ispc()
        [file,path] = uiputfile(file_suppoted,'Save',...
            filepath+"\"+fname);
    else
        [file,path] = uiputfile(file_suppoted,'Save',...
            filepath+"/"+fname);
    end
end

if isequal(file,0)
    status = SAVE_FAILED;
else
    filename = fullfile(path,file);
    [~,~,ext] = fileparts(filename);
    disp('data saving ...');
    s = size(I);
    I = reshape(I,[s(1:2),1,s(3:4)]);   % increase 4D to 5D
    switch ext
        case ".tif"
            if exist("metadata","var")
                % =====================
                bfCheckJavaMemory();
                bfCheckJavaPath();
                % must init path at first for calling java methods
                miniMetadata = createMinimalOMEXMLMetadata(I, dimOrder);
                fields = fieldnames(metadata);
                miniMetadata.setDatasetDescription(string(numel(fields)), 0);
                for k = 1:numel(fields)
                    data_store = string(fields{k}) + "=" + ...
                        string(metadata.(fields{k}));
                    miniMetadata.setDatasetDescription(data_store, k);
                end
                bfsave(I,filename,'Compression',compr,...
                    'BigTiff',true,'metadata',miniMetadata);
            else
                bfsave(I,filename,dimOrder,'Compression',compr,...
                    'BigTiff',true);
            end
        case ".mat"
            save(filename,"I","-mat");
        otherwise
    end
    
    status = SAVE_SUCCESS;
end
end

