classdef regpm < handle
    %REGPROJ This class defines registration project, which support <create
    %new project>, <open project>, <close project> and <save project>

    properties(Access = private, Constant)

    end

    % Register intecaction interface
    properties (GetAccess = ?Register, Dependent)
        Project             (1,1)       %
    end

    properties (GetAccess = public, SetAccess = immutable)
        OperationManager    (1,1)           % 1-by-1 regohm object
    end
    
    properties (GetAccess = private, SetAccess = immutable, Hidden)
        caller  (1,1)           % 1-by-1 Register object
    end

    properties (Access = private)
        project     (1,1)   regproj = regproj()
        opts        (1,1)   struct  = struct()
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

        function r = get.Project(this)
            r = this.project;       % shallow copy a handle
        end
    end

    methods (Access = public)
        % This function creates project on given location
        function [status, pf] = CreateNewProject(this, conf)
            arguments
                this
                conf    (1,1)   struct
            end

            % can't create new project if current project is active
            if ~isempty(this.project)
                status = -1;
                pf = "";
                return;
            else
                % create new empty regproj object
                this.project = regproj();

                % create project files by empty project
                pf = this.create_project_files(conf);
            end
        end

        function OpenProject(this, file)
            arguments
                this
                file    (1,1)   string  {mustBeFile}
            end

            
        end

        function CloseProject(this)

        end

        function SaveProject(this)

        end

        function r = FindAllProjects(this)
            
        end
    end

    methods (Access = private)
        % This function creates project folder on disk, which contains
        % files and folders
        function pf = create_project_files(this, conf)
            % Project Structure:
            % - <Data Folder>\*.dat     (project temporary data file)
            % - <Project Folder>
            %       |- [.Others]        (project misc. files)
            %       |- *.regproj        (main project file)
            %       |- *.log            (app log file)

            this.opts = conf;       % project creating options

            

        end
    end
end

