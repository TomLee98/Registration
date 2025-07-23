You can design your custom plugins for analysis in this folder, Register will 
check them when main app initializing and generate invoking interface in 
menu hierarchy "/Advanced/PluginsManager/". 
Note that Register only share some public variables after initialization. 
You need  to control the plugin to acquire data from Register by itself.

Note that the custom plugins need to follow the next 7 rules:

(1) A plugin must be a MATLAB executable file: *.m(function), *.p(function), *.mlapp(app) are supported.

(2) A plugin must be in a folder with following hierarchy (app as example):
    <my_plugin>
        |- my_plugin.mlapp
        |- my_plugin.xml
        |- <dependencies>
        |- <lang> (optional)

    where <dependencies> folder includes support files, my_plugin.xml is optional file, which includes
    information about my_plugin, such as "name", "version", "release date", "author(s)", etc. 
    You can package a plugin use .\plugins\PackageCustomPlugin.mlx for folder generation.
    Open PackageCustomPlugin.mlx, modify it, and run it by command window. (Do not use editor running!)
    Register will use "name" field in my_plugin.xml as plugin identifier under "PluginsManager" menu.
    And <lang> is an optional folder, you can put language mapping file .xml under here if multi-lanuage
    support is needed.

(3) A plugin must accept at least one input argument <caller> as 1-by-1 Register object, which must be 
    the first argument(for function) or the second argument(for app).

    (a) The function definition as:
        function [varargout] = my_plugin(caller, varargin)
            ...
        end
    
    (b) The app startup definition as:
        function startupFcn(app, caller)
            ...
        end

(4) A plugin can only visit the following shared properties (GetAccess = public) from Register.
    These variables are ('[',']' inner marked):

    % [Version] - Specifies the Register version
    % Read/Write Access - Read-Only
    % Accepted Values - string scalar, format as <main_ver>.<sub_ver>.<patch>, 
    %                                   which placeholder is oct nonnegtive integer
    % Default - N/A

(5) A plugin can only visit the following shared methods (Access = public) from Register.
    These methods are ('[',']' inner marked):

    % [GetImgInfo] - Return current image info
    % Usage: info = caller.GetImgInfo()
    % Input: None
    % Output: 1-by-12 image info table

    % [GetAppPosition] - Return main app figure  position(from left cornor), 
    %                    as [X, Y, W, H]
    % Usage: pos = caller.GetAppPosition()
    % Input: None
    % Output: 1-by-4 double array

    % [GetLangsInfo] - Return local machine region and languages mapping
    % Usage: [lang, mapping] = caller.GetLangsInfo()
    % Input: None
    % Output:
    %   - lang: 1-by-1 string, could be "zh_CN" or "en_US"
    %   - mapping: 1-by-1 struct, as language mapping (from Register basic dictionary)

    % [GetAppStyle] - Return app global style
    % Usage: stl = caller.GetAppStyle()
    % Input: None
    % Output:
    %   - stl: 1-by-1 string, could be "FOLLOW","LIGHT","DARK"

    % [GetAppTipsLevel] - Return app tips level
    % Usage: tl = caller.GetAppTipsLevel()
    % Input: None
    % Output:
    %   - tl: 1-by-1 string, could be "SIMPLE", "DETAIL", "NONE"

    % [GetAppMessageLevel] - Return app message 
    % Usage: ml = caller.GetAppMessageLevel()
    % Input: None
    % Output:
    %   - ml: 1-by-1 string, could be "ERROR", "WARNING", "INFO"

    % [PushInfo] - push plugin information into Register inner information collector
    % Usage: caller.PushInfo(subapp, obj, exception)
    % Input:
    %   - subapp: 1-by-1 string, current running plugin name
    %   - obj: 1-by-1 string, the name of function or callback which invokes [PushInfo]
    %   - exception: 1-by-1 string, the exception identifier or "None" if no exception
    % Output: None

    % [SendMessage] - send message to Register log that you want it is notified right now
    % Usage: caller.SendMessage(id, timeStamp)
    % Input:
    %   - id: 1-by-2 string, first item must be registered command (see 'Supplementary' 
              at the end of this file), second item could be parameter text
    %   - timeStamp: 1-by-1 logical, indicate if message with current time
    % Output: None

    % [SetProgressBar] - set the progress bar on Register panel
    % Usage: caller.SetProgressBar(x, iterinf)
    % Input:
    %   - x: 1-by-1 double in [0, 1], indicate the progress value
    %   - iterinf: 1-by-1 logical, if it is true and x is 0, progress bar will display
    %               an infinity loop animation
    % Output: None

(6) A plugin must take over its output, Register doesn't support outputs management.

(7) A plugin(*.mlapp) should implement function [Close], this will be
    called when Register close as parent.
    Note that these implement access permission need to be 'public', and
    visibility should not be hidden.

(8) A plugin could implement [RefreshPanelStyle], [RefreshPanelLanguage] 
    [RefreshMessageLevel] and [RefreshTooltips]. 
    Register could invoke these functions as needed.
    Note that these implement access permission need to be 'public', and
    visibility should not be hidden.

(9) A plugin shouldn't modify any Register public properties which don't 
    provide public implements, such as UI components.

======================================================================
You can follow the custom plugin design pipeline:

(1) Create a new plugin main app (or function),
(2) Put all dependencies into <dependencies> folder,
(3) Edit PackageCustomPlugin.mlx and RUN in command window.

Note that if you want to apply new plugin, please restart Register.

We design simple app view_max_activity.mlapp as an example.
Have fun :)




Supplementary:



