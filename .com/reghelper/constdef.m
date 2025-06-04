classdef constdef
    %CONSTDEF This class with some public constant defination
    
    properties (Access = public, Constant)
        GRAYSCALE_LIGHT = [0.94, 0.94, 0.94]
        GRAYSCALE_DARK =  [0.50, 0.50, 0.50]

        GRAYSCALE_AXES_LIGHT = [1.00, 1.00, 1.00]
        GRAYSCALE_AXES_DARK = [0.80, 0.80, 0.80]
    end

    properties (Access = {?mpimg, ?mpimgs}, Constant, Hidden)
        % code as: BUFFER_SIZE_MAX = 128
        BUFFER_KEY = "VsPkzE1Oc/bWrOEaYsjTaf8u6fpaRmKoGnbElf+MZX4="
    end
    
end