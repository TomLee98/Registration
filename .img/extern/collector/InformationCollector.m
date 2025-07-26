classdef InformationCollector < handle
    %INFORMATIONCOLLECTOR This class defines an information collector which
    %could generate encrypted operations log files and send to remote host
    %on LabNas if possible

    properties(Access = public, Constant)
        VER = "1.0.0"
    end

    properties(Access = private, Constant, Hidden)
        COMPONENTS_DELIMITER = " | "
    end
    
    properties(Access = private, Hidden)
        folder              (1,1)   struct  = struct("local","", "remote","")
        file                (1,1)   struct  = struct("local","", "remote","")
        is_online           (1,1)   logical
        is_running          (1,1)   logical = false
        scoder              (1,1)   
    end

    properties(Access = public, Dependent)
        Mode
        Status
    end
    
    methods
        function r = get.Mode(this)
            if this.is_online == true
                r = "Online";
            else
                r = "Offline";
            end
        end

        function r = get.Status(this)
            if this.is_running == true
                r = "Running";
            else
                r = "Stop";
            end
        end
    end

    methods (Access = public)
        function this = InformationCollector(run_mode)
            arguments
                run_mode        (1,1)   string  {mustBeMember(run_mode, ["ONLINE", "OFFLINE"])}
            end
            %INFORMATIONCOLLECTOR A constructor

            this.is_online = (run_mode == "ONLINE");
            this.scoder = aesobj();
        end

        function Start(this, folder)
            arguments
                this
                folder  (1,1)   string  {mustBeFolder}
            end

            if ~ this.is_running
                this.folder.local = folder;

                % create local file
                try
                    [status, this.file.local] = InformationCollector.make_envinfo_file(...
                        this.folder.local, this.COMPONENTS_DELIMITER);
                    if status ~= 0
                        warning("InformationCollector:invalidFolder", ...
                            "No valid folder for containing log file.");
                    end
                catch ME
                    rethrow(ME);
                end

                % make new log file on server
                if this.is_online == true
                    if ispc()
                        eval(this.scoder.decrypt(constdef.RIF_PC));
                    elseif isunix()
                        eval(this.scoder.decrypt(constdef.RIF_UNIX));
                    else
                        throw(MException("InformationCollector:unsupportPlatform", ...
                            "Only windows and unix are supported."));
                    end

                    % REMOTE_INFORMATION_FOLDER
                    this.folder.remote = REMOTE_INFORMATION_FOLDER;

                    try
                        [status, this.file.remote] = InformationCollector.make_envinfo_file(...
                            this.folder.remote, this.COMPONENTS_DELIMITER);
                        if status ~= 0
                            warning("InformationCollector:invalidFolder", ...
                                "No valid folder for containing log file.");
                        end
                    catch ME
                        rethrow(ME);
                    end
                end

                this.is_running = true;
            else
                warning("InformationCollector:invalidStart", ...
                    "An information collector is already running.");
            end
        end

        function Stop(this)
            this.is_running = false;
        end

        function Push(this, app, obj, exception)
            arguments
                this
                app         (1,1)   string
                obj         (1,1)   string
                exception   (1,1)   string  = "None"
            end

            if this.is_running == true
                % make text
                str = sprintf("%s%s%s%s%s%s%s", ...
                    string(datetime("now")), this.COMPONENTS_DELIMITER, ...
                    app, this.COMPONENTS_DELIMITER, ...
                    obj, this.COMPONENTS_DELIMITER, ...
                    exception);

                str_encrypted = this.scoder.encrypt(str);

                try
                    % write to local file
                    fid = fopen(this.file.local, "a+");
                    fprintf(fid, "%s\r\n", str_encrypted);
                    fclose(fid);

                    % write to remote file
                    if this.is_online == true
                        fid = fopen(this.file.remote, "a+");
                        fprintf(fid, "%s\r\n", str_encrypted);
                        fclose(fid);
                    end
                catch ME
                    rethrow(ME);
                end
            else
                warning("InformationCollector:invalidPush", ...
                    "Information Collector is not running.");
            end
        end
    end

    methods(Static, Hidden)
        function file_ = genfilename(folder_)
            rng('shuffle');     % make sure a random generation

            % generate file name code: 18 chars
            code_idx = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];

            % random suffix for avoiding same filename
            % and hidden the file information
            filename_ = char(code_idx(randperm(18)));

            if isstring(folder_), folder_ = folder_.char(); end

            file_ = [folder_, filesep, filename_, char(constdef.LOG_FILE_EXT)];
        end

        function [status, fout] = make_envinfo_file(folder, delimiter)
            status = 0;
            % try to access
            try
                % genearate fout name
                if isfolder(folder)
                    files_exists = struct2table(dir(folder));
                    files_exists = string(files_exists.name);
                    files_exists(~files_exists.contains(constdef.LOG_FILE_EXT)) = [];
                    while true
                        % generate new file name
                        fout = InformationCollector.genfilename(folder);
                        [~, fname, ext] = fileparts(fout);
                        fname = [fname, ext]; %#ok<AGROW>
                        if ~ismember(string(fname), files_exists)
                            break;
                        end
                    end
                else
                    status = -1;
                    return;
                end
            catch ME
                rethrow(ME);
            end

            % encrypt host name
            scoder = aesobj();
            [~, hname] = system('hostname');
            hname = strtrim(hname);
            hname_encryted = scoder.encrypt(hname);

            % open file and write basic information
            [~, hwinfo] = GetHardwareInfo();
            pe = pyenv;
            tbx_in = squeeze(string(struct2cell(ver)))';
            tbx_in = tbx_in(:, 1);  % take the first columns: name

            try
                % create file
                fid = fopen(fout, "a+");

                fprintf(fid, "########## Log File Automatically Generated by " + ...
                    "Information Collector ##########\r\n");
                fprintf(fid, "[Host Name] %s\r\n", hname_encryted);
                fprintf(fid, "[Creation Time] %s\r\n", string(datetime("now")));
                fprintf(fid, "[CPU Info] %s\r\n", ...
                    join(hwinfo.cpu, delimiter));
                fprintf(fid, "[GPU Info] %s\r\n", ...
                    join(hwinfo.gpu, delimiter));
                fprintf(fid, "[Memory Info] %s\r\n", ...
                    replace(hwinfo.mem, newline, delimiter));
                fprintf(fid, "[MATLAB Version] %s\r\n", version);
                fprintf(fid, "[Toolbox Installed] %s\r\n", ...
                    join(tbx_in, delimiter));
                fprintf(fid, "[Java Version] %s\r\n", version('-java'));
                fprintf(fid, "[Python Version] %s%s%s\r\n", ...
                    pe.Version, delimiter, pe.Home);
                fprintf(fid, "[Register Version] %s%s%s\r\n", Register.VER, ...
                    delimiter, Register.RELEASE_STATUS);
                fprintf(fid, "[InformationCollector Version] %s\r\n", this.VER);
                fprintf(fid, "\r\n########## Operations Queue Automatically Generated by " + ...
                    "Reg3D ##########\r\n");

                fclose(fid);
            catch ME
                rethrow(ME);
            end
        end
    end
end

