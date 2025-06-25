classdef constdef
    %CONSTDEF This class with some public constant definition
    
    properties (Access = public, Constant)
        REG3D_UPDATE_FILE_TITLE_ROW_SIZE = 2

        CAMERA_BACKGROUND = 100

        GRAYSCALE_LIGHT = [0.94, 0.94, 0.94]
        GRAYSCALE_DARK =  [0.50, 0.50, 0.50]

        GRAYSCALE_AXES_LIGHT = [1.00, 1.00, 1.00]
        GRAYSCALE_AXES_DARK = [0.80, 0.80, 0.80]

        COLOR_IMAGING_RED = [1.00, 0.00, 0.00; ...
                             0.00, 0.00, 0.00; ...
                             0.00, 0.00, 0.00]
        COLOR_IMAGING_GREEN = [0.00, 0.40, 0.00; ...
                               0.00, 1.00, 0.00; ...
                               0.00, 0.00, 0.00]
        COLOR_IMAGING_RED_GREEN = [0.71, 0.29, 0.00; ...
                                   0.00, 1.00, 0.00; ...
                                   0.00, 0.00, 0.00]
        COLOR_IMAGING_RED_CYAN = [1.00, 0.00, 0.00; ...
                                  0.00, 1.00, 0.00; ...
                                  0.00, 1.00, 0.00]
        COLOR_IMAGING_GREEN_MAGENTA = [1.00, 0.00, 0.00; ...
                                       0.00, 1.00, 0.00; ...
                                       1.00, 0.00, 0.00]

        UPDATE_TEXT_COLORMAP = dictionary("urgent",     "red", ...
                                          "new",        "blue", ...
                                          "adjust",     "green", ...
                                          "fix bugs",   "magenta")

        % operator definition
        OP_LOAD = "LOAD"
        OP_CROP = "CROP"
        OP_REGISTER = "REGISTER"
        OP_SEGMENT = "SEGMENT"

        % table as transmission optimal solution
        % 1 for load storage, 0 for real-time calculate recovery
        TRANSMISSION_STORAGE_TABLE = array2table([1, 1, 1, 0; ...
                                                  0, 1, 0, 0; ...
                                                  1, 1, 0, 0], ...
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
                                 "CachePolicy",           "PERFORMANCE", ...        % "PERFORMANCE"/"RESOURCES"/"BALANCE" 
                                 "BufferCleanTrigger",    "EXIT", ...               % "EXIT"/"RT"/"OFF"
                                 "Capacity",              64, ...                   % positive double scalar, 1 ~ 128
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
        % code as: HDD_BUFFER_SIZE_MAX = 128;
        HDD_BUFFER_KEY = "I6rlKEDWTeAmmC2/Yd311cVU1vbubYc8JQLKkaoMTmA="

        % code as: MEM_CACHE_SIZE_MAX = 64;
        MEM_CACHE_KEY = "Oav0Qo3zI0Qy+gxvtb6XscThqd7rEl6T/hwhDj3iySk="
    end
    
end