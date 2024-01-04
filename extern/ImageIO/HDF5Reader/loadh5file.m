function MovInfo = loadh5file( fName )
% loadh5file, load a H5 file and read data 
%   Input: file name (with full path)
%   Output: MovInfo, struct, image data with movie information

%get metadata 
[piezo,volumes,lasers,t,img_idx,res,acq_start_time]=loadh5metadata(fName);

% format the length data
piezo=piezo(1:length(img_idx));
volumes=volumes(1:length(img_idx));
lasers=lasers(1:length(img_idx),:);
imgTimes=zeros(max(volumes),1);

% use mean time of each volumn
for ii=1:max(volumes)
    imgTimes(ii)=mean(t(volumes==ii));
end

% color and z
which_lasers=find(any(lasers,1));
zstack_pos = 20*(piezo(volumes==1 & lasers(:,which_lasers(1))));
zdiff_size=diff(zstack_pos);

%initialize array
img1=initializeimgs(volumes,lasers,res,img_idx);

disp('Loading images...');
[img1] = getallvolumes(fName,lasers,piezo,img_idx,res,volumes,t,img1);
disp('Loading complete.');

%package data into a structure
MovInfo=struct;

MovInfo.fName = fName; % full path
MovInfo.sizeX = size(img1, 2);
MovInfo.sizeY = size(img1, 1);
MovInfo.sizeZ = length(zstack_pos);
MovInfo.sizeC = length(which_lasers);
MovInfo.sizeT = length(imgTimes);        %number of time, z and channel
MovInfo.tVol = imgTimes;                     %time of each volumn
MovInfo.img = img1;

end