function [] = ND2Close(FilePointer)

if ispc()
    libname = 'Nd2ReadSdk';
elseif isunix()
    libname = 'libNd2ReadSdk';
end

if exist('FilePointer','var')
    calllib(libname, 'Lim_FileClose', FilePointer);
end

end

