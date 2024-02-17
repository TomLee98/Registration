function data = loadims(file, opts, tspan, turbo)
%LOADIMS This function load ims files with tspan
% Input:
%   - file: 1-by-1 string, the data full file path
%   - opts: 1-by-12 table, images properties
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
%   - turbo: 1-by-1 logical, use external library for fast loading
% Output:
%   - info: struct with {opts, rt}, where opts is 1-by-12 metadata table,
%     rt is n-by-1 nonnegtive double time array
%
% where you can pass the time span as tspan for continous time cut

arguments
    file    (1,1)   string
    opts    (1,12)  table
    tspan   (1,:)   double = []
    turbo   (1,1)   logical = false
end

% generate HDF5 reader
GID = H5G.open(H5F.open(file), '/DataSet');
dataobj = DatasetReader(GID);

% load the image data
data = dataobj.GetData(tspan);

    function cOrder = getChannelOrder(cInfo)
        order = ["r","g","b"];
        cOrder = strings(1, numel(cInfo));
        for k = 1:numel(cInfo)
            [~, p] = max(cInfo(k).Color);
            cOrder(k) = order(p);
        end
    end

    function rt = getRelativeTime(timeStamp)
        rt = datetime(string(timeStamp),"Format","uuuu-MM-dd HH:mm:ss.SSS");
        rt = seconds(rt - rt(1));
    end
end

