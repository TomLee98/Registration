function [FilePointer, ImagePointer, ImageReadOut] = ND2Open(FileName)

if ispc()
    libname = 'Nd2ReadSdk';
elseif isunix()
    libname = 'libNd2ReadSdk';
end

if not(libisloaded(libname))
    warning('off', 'MATLAB:loadlibrary:cppoutput')
    [~, ~] = loadlibrary(libname, 'Nd2ReadSdk.h');
end

FileID = libpointer('voidPtr', [int8(FileName) 0]);
[FilePointer] = calllib(libname, 'Lim_FileOpenForReadUtf8', FileID);
Attibutes = calllib(libname, 'Lim_FileGetAttributes', FilePointer);
setdatatype(Attibutes, 'uint8Ptr', 500)
AttibutesValue = Attibutes.Value';
Attibuteslength = find(AttibutesValue == 0, 1);
AttibutesJson = char(AttibutesValue(1:Attibuteslength - 1));
AttibutesStru = jsondecode(AttibutesJson);

ImageStru.uiBitsPerComp = AttibutesStru.bitsPerComponentInMemory;
ImageStru.uiComponents = AttibutesStru.componentCount;
ImageStru.uiWidthBytes = AttibutesStru.widthBytes;
ImageStru.uiHeight = AttibutesStru.heightPx;
ImageStru.uiWidth = AttibutesStru.widthPx;

if ImageStru.uiWidthBytes==ImageStru.uiWidth*ImageStru.uiComponents*ImageStru.uiBitsPerComp/8
else
    warning('off','backtrace')
    warning('Image width is not fit the bytes of width. Reset image width.')
    warning('on','backtrace')
    ImageStru.uiWidth=ImageStru.uiWidthBytes/ImageStru.uiComponents/(ImageStru.uiBitsPerComp/8);
end

ImagePointer = libpointer('s_LIMPICTUREPtr', ImageStru);

calllib(libname, 'Lim_InitPicture', ImagePointer, ImageStru.uiWidth, ImageStru.uiHeight, ImageStru.uiBitsPerComp, ImageStru.uiComponents);

[~, ~, ImageReadOut] = calllib(libname, 'Lim_FileGetImageData', FilePointer, uint32(0), ImagePointer);

% explicit call member function
% Note that MATLAB R2025a Update-1 will throw error if implicit call this
% function
ImageReadOut.pImageData.setdatatype('uint16Ptr', ...
    ImageStru.uiWidth * ImageStru.uiHeight * ImageStru.uiComponents);
end
