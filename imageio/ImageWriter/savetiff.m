function [status, path] = savetiff(file, mov, ch, tspan, metadata, turbo)
%SAVETIFF This function save regmov as a tiff file (OME-TIFF or ImageJ-TIFF)
% input:
%   - file: 1-by-1 string, the exported file name
%   - mov: 1-by-1 regmov object
%   - ch: 1-by-c string, indicate the channel, could be "r", "g" or "b"
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - metadata: 1-by-1 struct, with field {ftype, operation}
%   - %   - turbo: 1-by-1 logical, use external library for fast saving
% output:
%   - status: the return status, SAVE_SUCCESS = 0; SAVE_FAILED = -1;
%   - path: the file saving folder
%
%   see also: bfsave

arguments
    file        (1,1)   string
    mov         (1,1)   regmov
    ch          (1,:)   string  {mustBeMember(ch, ["r","g","b"])}
    tspan       (1,2)   double  {mustBePositive, mustBeInteger} %#ok<INUSA> fully save
    metadata    (1,1)   struct
    turbo       (1,1)   logical = false
end

% memory size validation
mem = GetAvailableMemory();
if mem < (mov.Bytes.mem + mov.Bytes.map)
    warning("No enough memory for tiff saving.");
    status = -1;
    return;
end

[path, ~, ~] = fileparts(file);
file = char(file);

% extract data from regmov to memory
chs = mov.MetaData.cOrder;
img = mov.Movie(:,:,ch==chs,:,:);   % 5D array

if turbo == true
    pathToWriter = fileparts(mfilename('fullpath'));
    pathToSTPY = [pathToWriter, filesep, 'save_tiff.py'];
    if count(py.sys.path, pathToSTPY) == 0
        insert(py.sys.path,int32(0), pathToSTPY);
    end
    % append axes information in metadata
    metadata.axes = 'TZCYX';
    % append compression information
    metadata.compression = '';
    % generate a python dictionary
    mdata = py.dict(metadata);

    % permute I as imageJ hyperstack dimorder: 'TZCYXS', ordinary: 'XYCZT'
    img = permute(img, [5,4,3,1,2]);
    v = py.numpy.array(img).astype("uint16");
    fname = py.str(file);

    clearvars img
    try
        status = ...
            pyrunfile("save_tiff.py", "status", file=fname, vol=v, mdata=mdata);
    catch ME
        throwAsCaller(ME);
    end
else
    % using bio-formats package
    bfCheckJavaMemory();
    bfCheckJavaPath();
    % must init path at first for calling java methods
    miniMetadata = createMinimalOMEXMLMetadata(img, 'XYCZT');
    fields = fieldnames(metadata);
    miniMetadata.setDatasetDescription(string(numel(fields)), 0);
    for k = 1:numel(fields)
        data_store = string(fields{k}) + "=" + ...
            string(metadata.(fields{k}));
        miniMetadata.setDatasetDescription(data_store, k);
    end
    try
        bfsave(img, file, 'Compression','',...
            'BigTiff',true,'metadata',miniMetadata);
        status = 0;
    catch ME
        throwAsCaller(ME);
    end
end

end

