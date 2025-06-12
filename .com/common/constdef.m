classdef constdef
    %CONSTDEF This class with some public constant definition
    
    properties (Access = public, Constant)
        REG3D_UPDATE_FILE_TITLE_ROW_SIZE = 2

        CAMERA_BACKGROUND = 100

        GRAYSCALE_LIGHT = [0.94, 0.94, 0.94]
        GRAYSCALE_DARK =  [0.50, 0.50, 0.50]

        GRAYSCALE_AXES_LIGHT = [1.00, 1.00, 1.00]
        GRAYSCALE_AXES_DARK = [0.80, 0.80, 0.80]

        UPDATE_TEXT_COLORMAP = dictionary("urgent","red", "new","blue", ...
            "adjust","green", "fix bugs","magenta")

        % operator definition
        OP_LOAD = "LOAD"
        OP_CROP = "CROP"
        OP_REGISTER = "REGISTER"
        OP_SEGMENT = "SEGMENT"

        % table as transmission optimal solution
        % 1 for load storage, 0 for real-time calculate recovery
        TRANSMISSION_STORAGE_TABLE = array2table([1, 1, 1, 1; ...
                                                  0, 0, 1, 1; ...
                                                  1, 0, 1, 1], ...
            "VariableNames",["CROP",        "LOAD",     "REGISTER", "SEGMENT"], ...
            "RowNames",     ["PERFORMANCE", "RESOURCE", "BALANCE"]);

        NODE_FREE_TABLE = array2table(["KEEP",      "KEEP"; ...
                                       "DROP_ALL",  "KEEP"; ...
                                       "DROP_REG",  "KEEP"], ...
            "VariableNames",["OLDER",       "NEWER"], ...
            "RowNames",     ["PERFORMANCE", "RESOURCE", "BALANCE"]);

        PROFILE_DEFAULT = struct("Language",              "FOLLOW", ...             % "CN"/"US"/"FOLLOW"
                                 "SimpleStatistics",      "ON", ...                 % "ON"/"OFF"
                                 "AutoTemplate",          "ON", ...                 % "ON"/"OFF"
                                 "Style",                 "FOLLOW", ...             % "LIGHT"/"DARK"/"FOLLOW"
                                 "ProgressBarColor",      [0.30, 0.75, 0.93], ...   % 0 ~ 1, 1-by-3 array
                                 "PreferredOptimization", "PERFORMANCE", ...        % "PERFORMANCE"/"RESOURCES"/"BALANCE" 
                                 "CleanBufferTrigger",    "RT", ...                 % "EXIT"/"RT"/"OFF"
                                 "Capacity",              64, ...                   % positive double scalar, 1 ~ 128
                                 "BranchDepth",           "OFF", ...                % "LAST"/"LAST3"/"INFINITY"/"OFF"
                                 "NumProtectedCPU",       0, ...                    % nonnegative integer
                                 "NumProtectedGPU",       0, ...                    % nonnegative integer
                                 "MessageLevel",          "WARNING", ...            % "WARNING"/"INFO"
                                 "DataProtected",         "OFF", ...                % "ON"/"OFF"
                                 "PythonPath",            "", ...                   % string scalar as python executable path
                                 "AutoUpdate",            "RT", ...                 % "STARTUP"/"RT"/"EVERYDAY"/"OFF"
                                 "UpdateChannel",         "REL", ...                % "REL"/"PRE"/"ALL"
                                 "ExperimentalFeature",   "OFF")                    % "ON"/"OFF"
    end

    properties (Access = {?mpimg, ?mpimgs}, Constant, Hidden)
        % code as: BUFFER_SIZE_MAX = 128
        BUFFER_KEY = "VsPkzE1Oc/bWrOEaYsjTaf8u6fpaRmKoGnbElf+MZX4="
    end
    
end