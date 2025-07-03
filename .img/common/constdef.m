classdef constdef
    %CONSTDEF This class with some public constant definition
    
    properties (Access = public, Constant)
        %% environment definition
        SERVER_HOST_NAME = "silab"
        REGISTRATION_GROUP_NAME = "regusers"

        %% experiment definition
        CAMERA_BACKGROUND = 100

        %% updating definition
        REG3D_UPDATE_FILE_TITLE_ROW_SIZE = 2

        UPDATE_TEXT_COLORMAP = dictionary("urgent",     "red", ...
                                          "new",        "blue", ...
                                          "adjust",     "green", ...
                                          "fix bugs",   "magenta")

        %% appearance definition
        GRAYSCALE_LIGHT = [0.94, 0.94, 0.94]
        GRAYSCALE_DARK =  [0.50, 0.50, 0.50]

        GRAYSCALE_AXES_LIGHT = [1.00, 1.00, 1.00]
        GRAYSCALE_AXES_DARK = [0.80, 0.80, 0.80]

        %% colormap definition
        COLOR_PERCENTAGE_SAFE       = [0.25,0.75,1.00]
        COLOR_PERCENTAGE_URGENT     = [1.00, 0.25, 0.25]

        COLOR_IMAGING_RED           = [1.00, 0.00, 0.00; ...
                                       0.00, 0.00, 0.00; ...
                                       0.00, 0.00, 0.00]

        COLOR_IMAGING_GREEN         = [0.00, 0.40, 0.00; ...
                                       0.00, 1.00, 0.00; ...
                                       0.00, 0.00, 0.00]

        COLOR_IMAGING_RED_GREEN     = [0.71, 0.29, 0.00; ...
                                       0.00, 1.00, 0.00; ...
                                       0.00, 0.00, 0.00]

        COLOR_IMAGING_RED_CYAN      = [1.00, 0.00, 0.00; ...
                                       0.00, 1.00, 0.00; ...
                                       0.00, 1.00, 0.00]

        COLOR_IMAGING_GREEN_MAGENTA = [1.00, 0.00, 0.00; ...
                                       0.00, 1.00, 0.00; ...
                                       1.00, 0.00, 0.00]

        COLOR_IMAGIMG_GRAY          = [1.00, 0.00, 0.00; ...
                                       0.00, 1.00, 0.00; ...
                                       0.00, 0.00, 1.00]

        %% merge mode definition
        SOURCE_MERGE_NORMAL = "NORMAL"
        SOURCE_MERGE_PREPROCESS = "PRE_PROCESS"
        SOURCE_MERGE_FUSED = "FUSED"

        %% restored operation definition
        OP_LOAD = "LOAD"
        OP_CROP = "CROP"
        OP_REGISTER = "REGISTER"
        OP_SEGMENT = "SEGMENT"

        %% cache policy definition
        % table as transmission optimal solution
        % 1 for load storage, 0 for real-time calculate recovery
        TRANSMISSION_STORAGE_TABLE = array2table([1, 1, 1, 0; ...
                                                  0, 1, 0, 0; ...
                                                  1, 1, 0, 0], ...
            "VariableNames",["CROP",        "LOAD",     "REGISTER", "SEGMENT"], ...
            "RowNames",     ["PERFORMANCE", "RESOURCE", "BALANCE"]);

        %% app profile definition
        PROFILE_DEFAULT = struct("Language",              "FOLLOW", ...             % "CN"/"US"/"FOLLOW"
                                 "SimpleStatistics",      "ON", ...                 % "ON"/"OFF"
                                 "AutoTemplate",          "ON", ...                 % "ON"/"OFF"
                                 "PreprocessEvaluation",  "ON", ...                 % "ON"/"OFF"
                                 "Style",                 "FOLLOW", ...             % "LIGHT"/"DARK"/"FOLLOW"
                                 "StorageViewStyle",      "MODERN", ...             % "CLASSICAL"/"MODERN"
                                 "ProgressBarColor",      [0.30, 0.75, 0.93], ...   % 0 ~ 1, 1-by-3 array
                                 "CacheLocation",         "AUTO", ...               % "AUTO"/"CUSTOMIZED"
                                 "CachePolicy",           "PERFORMANCE", ...        % "PERFORMANCE"/"RESOURCES"/"BALANCE" 
                                 "CacheCleanTrigger",     "EXIT", ...               % "EXIT"/"OFF"
                                 "MemoryCapacity",        64, ...                   % positive double scalar, 1 ~ 96
                                 "HardDriveCapacity",     128, ...                   % positive double scalar, 1 ~ 256
                                 "NumProtectedCPU",       0, ...                    % nonnegative integer
                                 "NumProtectedGPU",       0, ...                    % nonnegative integer
                                 "MessageLevel",          "WARNING", ...            % "WARNING"/"INFO"
                                 "DataProtected",         "OFF", ...                % "ON"/"OFF"
                                 "PythonPath",            "", ...                   % string scalar as python executable path
                                 "AutoUpdate",            "RT", ...                 % "STARTUP"/"RT"/"EVERYDAY"/"OFF"
                                 "UpdateChannel",         "REL", ...                % "REL"/"PRE"/"ALL"
                                 "ExperimentalFeature",   "OFF")                    % "ON"/"OFF"
    end

    properties (Access = {?mpimg, ?mpimgs, ?regohm}, Constant, Hidden)
        % code as: HDD_BUFFER_SIZE_MAX = 256;
        HDD_BUFFER_KEY = "I6rlKEDWTeAmmC2/Yd311bJ+f9h+K+9CfCYICfRZKcY="

        % code as: MEM_CACHE_SIZE_MAX = 96;
        MEM_CACHE_KEY = "Oav0Qo3zI0Qy+gxvtb6XscbuMv0sbTJiUvNQY9m9kRo="
    end
    
end