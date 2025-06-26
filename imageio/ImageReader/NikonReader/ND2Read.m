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
        % explicit call member function
        % Note that MATLAB R2025a Update-1 will throw error if implicit call this
        % function
        
        Image = reshape(ImageReadOut.pImageData, [ImageReadOut.uiComponents, ImageReadOut.uiWidth * ImageReadOut.uiHeight]);

        for j = 1:ImageReadOut.uiComponents
            ImageStack{j}(:, :, i) = reshape(Image(j, :), [ImageReadOut.uiWidth, ImageReadOut.uiHeight])';
        end

    end

end

end
