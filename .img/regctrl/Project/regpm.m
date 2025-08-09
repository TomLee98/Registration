classdef regpm < handle
    %REGPROJ This class defines registration project, which support <create
    %new project>, <open project>, <close project> and <save project>

    % Register intecaction interface
    properties (Access = ?Register, Dependent)
        IsProjectActive     (1,1)       % get/___, 1-by-1 logical, indicate if project on hard-drive
        Project             (1,1)       % get/___, 1-by-1 regproj object
        ProjectName         (1,1)       % get/___, 1-by-1 string
    end

    properties (GetAccess = public, SetAccess = immutable)
        OperationManager    (1,1)           % 1-by-1 regohm object
    end
    
    properties (GetAccess = private, SetAccess = immutable, Hidden)
        caller  (1,1)           % 1-by-1 Register object
    end

    properties (Access = private)
        project     (1,1)   regproj = regproj.empty()
        proj_opts   (1,1)   struct  = struct()
        active_flag (1,1)   logical = false
    end
    
    methods
        %% Constructor
        function this = regpm(caller, opts)
            %REGPROJ A Constructor
            arguments
                caller  (1,1)   Register
                opts    (1,1)   struct
            end

            % set caller
            this.caller = caller;

            % initialize ohmgr
            this.OperationManager = regohm(caller.GetUIFigure(), ...
                                           caller.Tree_Operations, ...
                                           caller.TextArea_OperationSnapshot, ...
                                           caller.StorageManager, ...
                                           caller.ContextMenu_NodeOperation, ...
                                           opts.MemoryCapacity, ...
                                           opts.HardDriveCapacity, ...
                                           opts.CachePolicy);
        end
        
        %% Getter/Setter
        function r = get.IsProjectActive(this)
            r = this.active_flag;
        end

        function r = get.Project(this)
            r = this.project;       % shallow copy a handle
        end

        function r = get.ProjectName(this)
            r = this.project.proj_fname;
        end
    end

    methods (Access = public)
        % This function creates project on given location
        function [status, pf, dst] = CreateNewProject(this, prof, conf)
            arguments
                this
                prof    (1,1)   struct
                conf    (1,1)   struct
            end

            % can't create new project if current project is active
            if this.active_flag == true
                status = -1;
                pf = "";
                dst = "";
                return;
            else
                % create new empty regproj object
                this.project = regproj.empty();

                % create project files by empty project
                [pf, dst] = this.create_project_files(prof, conf);

                % modify the cache policy
                this.OperationManager.CachePolicy = conf.CachePolicy;

                this.project.IsSaved = true;
                this.active_flag = true;
                status = 0;
            end
        end

        function status = OpenProject(this, pfile)
            arguments
                this
                pfile    (1,1)   string  {mustBeFile}
            end

            % can't open project if current project is active
            if this.active_flag == true
                status = -1;
                return;
            else
                % create new empty regproj object
                this.project = regproj.empty();

                % create project files by empty project
                this.open_exist_project(pfile);

                % this.project.IsSaved = true;
                this.active_flag = true;
                status = 0;
            end
        end

        function CloseProject(this)
            % clear temporary files (which could be auto generated if project is reopened)
            % ...
            
            % delete project (but keep the files)
            try
                % clear operation history and manager
                this.OperationManager.free();

                % clear project
                occupied = false;
                save(this.ProjectName, "occupied", '-mat', '-append');
                
                delete(this.project);
                this.project = regproj.empty();   % reset as an empty project

                this.active_flag = false;
            catch ME
                rethrow(ME);
            end
        end

        function SaveProject(this)
            if isfile(this.project.proj_fname)
                this.project.IsSaved = true;
                DS_ = this.project.Seed;    % get data struct
                try
                    % save project data
                    save(this.project.proj_fname, "DS_", '-mat', '-v7.3', '-nocompression');

                    % save operation history
                    this.OperationManager.Save(this.project.proj_fname);
                catch ME
                    rethrow(ME);
                end
            end
        end

        function ResetProject(this)
            % This function resets project variables as initialized
            if ~isempty(this.project)
                %
                this.project.Reset();

                % keep project option
            end
        end
    end

    methods (Access = private)
        % This function creates project folder on disk, which contains
        % files and folders
        function [pfile, dst] = create_project_files(this, prof, conf)
            % Project Structure:
            % - <Data Folder>\*.dat     (project temporary data file)
            % - <Project Folder>
            %       |- [.misc]          (project misc. files)
            %       |- *.regproj        (main project file)
            %       |- *.log            (app log file, create by Logger)

            this.proj_opts = conf;       % project creating options
            pfolder = this.proj_opts.ProjectFolder;
            pfile = pfolder + filesep + this.proj_opts.ProjectName + constdef.PROJECT_FILE_EXT;

            %% create project folder
            try
                % remake dir
                if isfolder(pfolder), rmdir(pfolder, 's'); end
                mkdir(pfolder);
            catch ME
                rethrow(ME);
            end

            %% put project files
            % 1. Create subfolder
            % this folder contains mask file and CaImAn initial value file
            optree_folder = pfolder + filesep + constdef.OPTREE_FOLDER_NAME;
            try
                mkdir(optree_folder);
            catch ME
                rethrow(ME);
            end

            % 2. Make project file: *.regproj
            try
                % save data field
                DS_ = this.project.Seed;    % get data struct
                occupied = true;
                save(pfile, "DS_", "occupied", '-mat', '-v7.3', '-nocompression');
                % save operation field
                this.OperationManager.Save(pfile);
            catch ME
                rethrow(ME);
            end

            %% create data folder
            dst = prof.CacheRootFolder + filesep + this.proj_opts.CacheFolder;
            try
                % remake dir
                if isfolder(dst), rmdir(dst, 's'); end
                mkdir(dst);
            catch ME
                rethrow(ME);
            end
        end

        function open_exist_project(this, pfile)
            % construct all project from exist file
            try
                % load project data from file
                load(pfile, '-mat', "DS_");
                this.project.Restore(DS_);

                % load operation history from file
                this.OperationManager.Load(pfile);

                occupied = true;
                save(pfile, "occupied", '-mat', '-append');
            catch ME
                rethrow(ME);
            end
        end
    end
end

