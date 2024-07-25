function metadata = LoadMetaData(time)
%LOADMETADATA This function load the experiment trial to compare with the 
% time vector and output the string vector which has the same length as
% time
% Input:
%   - time: double, n*1 time vector
% Output:
%   - metadata: 1*n cell

% Version: 1.0.0
%   *** Support bad valve skipped

arguments
    time (:,1) double;  % time vector
end

tpn = numel(time);
metadata = cell(1, tpn);

% load the metadata file
pathToMFile = fileparts(mfilename('fullpath'));
[file, path] = uigetfile({'*.xlsx','Excel 文件 (*.xlsx)';...
    '*.xls','Excel 2003 文件 (*.xls)';...
    '*.csv','逗号分隔符文件 (*.csv)'},...
    "打开元数据(气味-时间标记)",...
    pathToMFile+"\metadata.xlsx");
if path ~= 0
    file = fullfile(path, file);
    data = readtable(file);
    % remove the possible none rows
    data = modifysdf(data);

    odors = string(table2array(data(3:end,1)));
    odor_seq = reshape(table2array(data(1:2,3:end)),[],1);
    odor_seq = cumsum([0;odor_seq]);
    odor_mixtable = logical(table2array(data(3:end,3:end)));
    for k = 1:tpn
        p = find(time(k)<odor_seq,1,"first");
        if ~isempty(p)
            switch mod(p,2)
                case 0
                    % drop in water
                    metadata{k} = "water";
                case 1
                    % drop in odors
                    metadata{k} = odors(odor_mixtable(:,(p-1)/2));
            end
        else
            % meta data lost
            metadata{k} = "none";
        end
    end
else
    metadata = [];
end

end

