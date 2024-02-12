function [piezo,volumes,lasers,t,imgIdx,res,startTime] = loadh5metadata(fName)
% loadh5imgproperties, load the metadata of a h5 file
%   Input: file name (with full path)
%   Output: metadata

    imgIdx=double(h5read(fName,'/img_idx'));
    dAQ=h5read(fName,'/DAQ');
    res=h5readatt(fName,'/','Resolution');
    t=double(h5read(fName,'/t'))/1e9;
    piezo=dAQ(imgIdx,1);
   
    lasers=dAQ(imgIdx,2:end-1);
    volumes=dAQ(imgIdx,end);
    startTimeStr = h5readatt(fName,'/','Start Time');
    startTime = datetime(startTimeStr,'InputFormat','yyyyMMdd HH:mm:ss');

%     start_time=datetime(start_time_str{1},'InputFormat','yyyyMMdd HH:mm:ss');
