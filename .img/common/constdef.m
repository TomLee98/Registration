classdef constdef
    %CONSTDEF This class with some public constant definition
    
    properties (Access = public, Constant)
        %% basic definition
        KB2B = 1024
        MB2B = 1048576
        GB2B = 1073741824

        %% environment definition
        SERVER_HOST_NAME = "silab"
        REGISTRATION_GROUP_NAME = "regusers"
        MEX_SETUP = "MATLAB_MEX_SETUP"
        REGPROGJ_REGISTRY_FLAG = "Reg3Dfile"
        PROFILE_FOLDER = ['Reg3D', filesep]
        APP_PROFILE_FILE_NAME = ".profile.xml"
        PROJECT_PROFILE_FILE_NAME = ".project.xml"
        PROJECT_NAME_DEFAULT = "Untitled"
        PROJECT_FOLDER_NAME_DEFAULT = "Reg3D Projects"
        OPTREE_FOLDER_NAME = ".misc"
        SUBAPP_CLOSE_IMPLEMENT_NAME = "Close"
        SUBAPP_UPDATE_STYLE_IMPLEMENT_NAME = "RefreshPanelStyle"
        SUBAPP_UPDATE_LANGUAGE_IMPLEMENT_NAME = "RefreshPanelLanguage"
        SUBAPP_UPDATE_MESSAGE_IMPLEMENT_NAME = "RefreshMessageLevel"
        SUBAPP_UPDATE_TOOLTIP_IMPLEMENT_NAME = "RefreshTooltips"

        %% file extension definition
        RAW_FILE_EXT = ".dat"
        LOG_FILE_EXT = ".log"
        MASK_FILE_EXT = ".mask"
        CMMASK_FILE_EXT = ".mat"
        PARALLEL_PROFILE_EXT = ".mlsettings"
        PROJECT_FILE_EXT = ".regproj"
        REGMOV_FILE_EXT = ".rmv"
        TIME_FILE_EXT = ".tim"
        HLOG_FILE_EXT = ".txt"
        APPCONF_FILE_EXT = ".xml"
        LANGUAGE_FILE_EXT = ".xml"

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

        %% tooltips level definition
        TOOLTIP_SIMPLE  = "SIMPLE"
        TOOLTIP_DETAIL  = "DETAIL"
        TOOLTIP_NONE    = "NONE"

        %% message level definition
        MESSAGE_ERROR   = "ERROR"
        MESSAGE_WARNING = "WARNING"
        MESSAGE_INFO    = "INFO"

        %% restored operation definition
        OP_LOAD = "LOAD"
        OP_CROP = "CROP"
        OP_REGISTER = "REGISTER"
        OP_SEGMENT = "SEGMENT"

        %% cache policy definition
        % table as transmission optimal solution
        % 1 for load storage, 0 for real-time calculate recovery
        TRANSMISSION_STORAGE_TABLE = array2table([1, 0, 1, 0; ...
                                                  0, 0, 0, 0; ...
                                                  1, 0, 0, 0], ...
            "VariableNames",["CROP",        "LOAD",     "REGISTER", "SEGMENT"], ...
            "RowNames",     ["PERFORMANCE", "RESOURCE", "BALANCE"]);

        %% app profile definition
        PROFILE_DEFAULT = struct("Language",              "FOLLOW", ...         % "CN"/"US"/"FOLLOW"
                                 "Style",                 "FOLLOW", ...         % "LIGHT"/"DARK"/"FOLLOW"
                                 "StorageViewStyle",      "MODERN", ...         % "CLASSICAL"/"MODERN"
                                 "ProgressBarColor",      [0.30, 0.75, 0.93],...% 0 ~ 1, 1-by-3 array
                                 "SimpleStatistics",      "ON", ...             % "ON"/"OFF"
                                 "AutoTemplate",          "ON", ...             % "ON"/"OFF"
                                 "PreprocessEvaluation",  "ON", ...             % "ON"/"OFF"
                                 "CacheRootFolder",       "", ...               % cache files root folder, under <user> data folder as usual
                                 "AutoSave",              "ON", ...             % "ON"/"OFF", if project is auto save
                                 "AutoSaveInterval",      10, ...               % positive integer scalar, 1 ~ 120, unit as minute
                                 "ImportSourcePreset",    "HOME", ...           % "LAST"/"MODE"/"HOME", indicate default import source
                                 "MemoryCapacity",        64, ...               % positive double scalar, 1 ~ 96
                                 "HardDriveCapacity",     128, ...              % positive double scalar, 1 ~ 256
                                 "NumProtectedCPU",       0, ...                % nonnegative integer
                                 "NumProtectedGPU",       0, ...                % nonnegative integer
                                 "MessageLevel",          "WARNING", ...        % "ERROR"/"WARNING"/"INFO"
                                 "TooltipsLevel",         "SIMPLE", ...         % "SIMPLE"/"DETAIL"/"NONE"
                                 "InformationCollection", "ON", ...             % "ON"/"OFF"
                                 "PythonPath",            "", ...               % string scalar as python executable path
                                 "DeveloperMode",         "OFF", ...            % "ON"/"OFF"
                                 "AutoUpdate",            "RT", ...             % "STARTUP"/"RT"/"EVERYDAY"/"OFF"
                                 "UpdateChannel",         "REL", ...            % "REL"/"PRE"/"ALL"
                                 "ExperimentalFeature",   "OFF")                % "ON"/"OFF"

        %% project profile defination (global)
        PROJECT_PROFILE_DEFAULT = struct("ProjectFolders",  "", ...             % string array as project root folder
                                         "SourceFolders",   "")                 % string array as import source folders

        %% project configuration (local)
        PROJECT_CONFIG_DEFAULT = struct("CachePolicy",      "PERFORMANCE", ...  % "PERFORMANCE"/"RESOURCES"/"BALANCE"
                                        "CacheFolder",      "", ...             % string scalar indicate data folder name
                                        "DataProtected",    "OFF", ...          % "ON"/"OFF"
                                        "ProjectFolder",    "", ...             % string scalar indicate current project folder
                                        "ProjectName",      "Untitled", ...     % "Untitled" as constant, can not be changed
                                        "Template",         "NONE", ...         % "NONE", "STDREG", "MEDIA", "DEV"
                                        "WaveLength",       [694, 525, 440])    % "r": 694nm, "g": 525nm, "b": 440nm as default
    end

    properties (Access = {?mpimg, ?mpimgs, ?regohm}, Constant, Hidden)
        % code as: HDD_BUFFER_SIZE_MAX = 256;
        HDD_BUFFER_KEY = "I6rlKEDWTeAmmC2/Yd311bJ+f9h+K+9CfCYICfRZKcY="

        % code as: MEM_CACHE_SIZE_MAX = 96;
        MEM_CACHE_KEY = "Oav0Qo3zI0Qy+gxvtb6XscbuMv0sbTJiUvNQY9m9kRo="
    end

    properties (Access = ?InformationCollector, Constant, Hidden)
        % code as: REMOTE_INFORMATION_FOLDER = '\\LabNas1\group_SharedFolder\Code\Reg3D App\.reports';
        RIF_PC = "S9DNLUrQeOcvdj+W5tSBW6xQyL9L1QKd9mAWLONqiI23uKl4tcsBL0uCWHNF5ptTFtK9btWWDKFuEtiRe5i7qNfbcPCY3bXLsgvyhv8ZVI4bzdoMnNcufVnJz1V5gAeV";

        % code as: REMOTE_INFORMATION_FOLDER = '/data/Share/Code/Reg3D App/.reports';
        RIF_UNIX = "S9DNLUrQeOcvdj+W5tSBW7zjwXUSstMziN530gTIkDTqkWsYR7BicaYvZy2OU17cMwIIiGj1ZBR7abrE+ClTdZhfo0tRX3ztnkwQqHarlSA=";
    end
    
end