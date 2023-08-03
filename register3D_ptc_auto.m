function [signal, tform, marker] = register3D_ptc_auto(varargin)
%REGISTER3D_PTC This function uses point clouds and registration chain do
%long term 3d registration
%
%   signal = register3D()
%   [signal,tform,marker] = register3D()
%   [signal,tform,marker] = register3D(filename)
%   [signal,tform,marker] = register3D(movinfo,filename)
%   [signal,tform,marker] = register3D(___,Name,Value)
%
%   input:
%   - (optional) movinfo:  the structure of movie data(4/5-D) and
%                          information structure with field {'mov','opts'}
%   - (optional) filename: the name of file need to be aligned, only *.ims,
%                          *.nd2 and *.tiff is valid, image dimension: x*y*(c)*z*t
%   - varargin:
%       - RefVol: the fixed volume indicator, struct with {'G','L'}, 'G' for
%                 'global' and 'L' for 'local', each element is 1*2 string,
%                 for [method, range]. method indicates Ref calculate method,
%                 range indicates the sampling frames.
%                 Example: RefVol = struct("G",["mean","1,2,3"], "L",["mean","1,2,3"]);
%       - RegFrames: the frames series need to be registered
%       - KeyFrames: the key frames series, as 4D tensor
%       - IntensityThreshold: the lower intensity threshold for
%                           pre-foreground extraction, 1~100, 97 as default
%       - ScaleThreshold: the scale threshold for pre-foreground
%                        extraction, omit little objects for robust registration
%       - DS_PointCloud: the downsampling options for point cloud
%                      registration, could be "GA","RAND","NU", "GA" as default
%       - DS_PointCloud_Param: the downsampling parameter value
%       - OutlierRatio: the outlier ratio for more anti-noise in points 
%                       cloud registration, 0~1, 0.1 as default
%       - MaxIterN_PointCloud: the maximum iteration number when points 
%                              cloud registration is running, which is one 
%                              of the stop condition, 50 as default
%       - ErrorLimit: the error limitation when points cloud is running,
%                     which is one of the stop condition
%       - InitStep: the initial step when optimizal two volumes loss, it is
%                   initial radius(default: 6.25e-3) when using "multimodal", or
%                   maximum step length(default: 6.25e-2) when using "monomodal"
%       - MinStep: the minimum iteration step when optimizal two volumes
%                  loss, which meature the convergance
%       - DS_Voxel: the downsampling options for speed up and more robust, could
%                   be "auto", "2X2", "3X3", "4X4", "auto" as default
%       - MaxIterN_Voxel: the maximum iteration number when alignment, 100 (default)
%       - IterCoeff: the iteration coefficient for loop, if Modal is monomodal,
%                    then the range of IterCoeff is (0,1), 0.5 as default;
%                    and if Modal is multimodal, then the range of
%                    IterCoeff is (0,+inf), 1.05 as default
%       - AutoContrast: the flag for auto contrast on marker channel, for
%                       better local registration, true(default)/false
%       - CompAcc: the compensate accuracy, the imhistmatchn binning counts
%       - AFS: the accumulative field smoothing, the AFS bigger and the
%              displacement is more smooth, usual 0.5-3, 1 as default
%       - Chla: the channel you want to align, this section is enabled
%                only if multi-channel ["r"(default),"g","b"]
%       - Chls: the channel you want to extract signal ["r","g"(default),"b"]
%       - LRTDS: the local registration transformation downsampling level,
%                integer number for how many pixel cause a sampling point,
%                1(no downsample) as default
%   output:
%   - tform: 1st-col: the global transform object, affine2d or affine3d
%            2nd-col: the local displacement field, double
%   - signal: the signal channel has been aligned, 4-D uint16
%   - marker: the marker channel has been aligned, 4-D uint16
%
%   see also: imregdemons, imregconfig, imregtform, imref3d, pcdownsample,
%             pcregistercpd, imwarp


end

