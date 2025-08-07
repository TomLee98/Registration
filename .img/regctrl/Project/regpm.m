classdef regpm < handle
    %REGPROJ This class defines registration project, which support <create
    %new project>, <open project>, <close project> and <save project>

    % Register intecaction interface
    properties (Access = ?Register, Dependent)
        Project             (1,1)       %
        ProjectName         (1,1)
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
                                           opts.HardDriveCapacity);
        end
        
        %% Getter/Setter
        function r = get.Project(this)
            r = this.project;       % shallow copy a handle
        end

        function r = get.ProjectName(this)
            if ~isempty(this.project)
                r = this.project.proj_fname;
            else
                r = "";
            end
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
            if ~isempty(this.project)
                status = -1;
                pf = "";
                dst = "";
                return;
            else
                % create new empty regproj object
                this.project = regproj.empty();

                % create project files by empty project
                [pf, dst] = this.create_project_files(prof, conf);

                this.project.IsSaved = true;

                status = 0;
            end
        end

        function [status, pf] = OpenProject(this, file)
            arguments
                this
                file    (1,1)   string  {mustBeFile}
            end

            if ~isempty(this.project)
                status = -1;
                pf = "";
                return;
            else
                % create new empty regproj object
                this.project = regproj.empty();

                % create project files by empty project
                pf = this.open_exist_project(file);

                status = 0;
            end
        end

        function CloseProject(this)
            % clear temporary files (which could be auto generated if project is reopened)
            % ...
            
            % delete project (but keep the files)
            if ~isempty(this.project)
                try
                    % clear operation history and manager
                    this.OperationManager.free();

                    % clear project
                    occupied = false;
                    save(this.ProjectName, "occupied", '-mat', '-append');
                    delete(this.project);
                    this.project = regproj.empty();   % reset as an empty project
                catch ME
                    rethrow(ME);
                end
            end
        end

        function SaveProject(this)
            if ~isempty(this.project)
                pfile = this.project.proj_fname;
                this.project.IsSaved = true;
                DS_ = this.project.Seed;    % get data struct
                try
                    % save project data
                    save(pfile, "DS_", '-mat', '-v7.3', '-nocompression');

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
                DS_ = this.project.Seed;    % get data struct
                occupied = true;
                save(pfile, "DS_", "occupied", '-mat', '-v7.3', '-nocompression');
            catch ME
                rethrow(ME);
            end

            %% create data folder
            dst = prof.DataRootFolder + filesep + this.proj_opts.DataFolder;
            try
                % remake dir
                if isfolder(dst), rmdir(dst, 's'); end
                mkdir(dst);
            catch ME
                rethrow(ME);
            end
        end

        function pfile = open_exist_project(this, file)
            % construct all project from exist file
            try
                % load project data from file
                load(file, '-mat', "DS_");
                this.project.Restore(DS_);

                % load operation history from file
                this.OperationManager.Load(file);

                occupied = true;
                save(file, "occupied", '-mat', '-append');
            catch ME
                rethrow(ME);
            end

            pfile = this.project.proj_fname;
        end
    end
end

