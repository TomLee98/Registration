function [ImageStack] = ND2Read(FilePointer, ImagePointer, ImageReadOut, Num)
if ispc()
    libname = 'Nd2ReadSdk';
elseif isunix()
    libname = 'libNd2ReadSdk';
end

if ImageReadOut.uiComponents == 1
    ImageStack = zeros([ImageReadOut.uiHeight, ImageReadOut.uiWidth, size(Num, 2)],'uint16');

    for i = 1:size(Num, 2)
        [~, ~, ImageReadOut] = calllib(libname, 'Lim_FileGetImageData', FilePointer, uint32(Num(i) - 1), ImagePointer);
        Image = reshape(ImageReadOut.pImageData, [ImageReadOut.uiWidth, ImageReadOut.uiHeight]);
        ImageStack(:, :, i) = Image';
    end

else
    ImageStack = cell([1, ImageReadOut.uiComponents]);

    for i = 1:size(Num, 2)
        [~, ~, ImageReadOut] = calllib(libname, 'Lim_FileGetImageData', FilePointer, uint32(Num(i) - 1), ImagePointer);
        
        % Note that if MATLAB version <= R2024b, ImageReadOut.pImageData
        % will be uint16 array, but if not, it will be libpointer
        if isMATLABReleaseOlderThan("R2025a")
            % call MATLAB reshape (overload for uint16 array)
            Image = reshape(ImageReadOut.pImageData, ...
                [ImageReadOut.uiComponents, ImageReadOut.uiWidth * ImageReadOut.uiHeight]);
        else
            % must explicit call member function!!!
            ImageReadOut.pImageData.setdatatype('uint16Ptr', ...
                ImageReadOut.uiWidth * ImageReadOut.uiHeight * ImageReadOut.uiComponents);
            ImageReadOut.pImageData.reshape(ImageReadOut.uiComponents, ...
                ImageReadOut.uiWidth * ImageReadOut.uiHeight);
            Image = ImageReadOut.pImageData.Value;
        end

        

        for j = 1:ImageReadOut.uiComponents
            ImageStack{j}(:, :, i) = reshape(Image(j, :), [ImageReadOut.uiWidth, ImageReadOut.uiHeight])';
        end

    end

end

end
