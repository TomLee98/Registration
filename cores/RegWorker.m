classdef RegWorker < handle
    %REGWORKER This class is regworker object defination
    % The RegWorker object can select registration algorithm and use the
    % TaskManager for multi-user calculation scheduling

    % user -> Run -> TaskManager -> align_one_vf

    properties(Hidden, Constant)
        WORKER_STATE = enumeration('matlab.lang.OnOffSwitchState')

        STATUS_SUCCESS = 0
        STATUS_FAILED = -1
    end

    properties(Access = private)
        regopt              % registration options
        volopt              % image volumes options

        mov_aligned         % the aligned mov
        mov_aligned_zproj   % the aligned mov z-projection
        mov_aligned_tproj   % the aligned mov t-projection
        tform               % the transform, n-by-3 cell
    end

    properties(Access = private, Hidden)
        caller      % caller, must be Register object
        taskmgr     % TaskManager Object
        state       % 1-by-1 OnOffSwitchState enum, can be "on" or "off" 
    end
    
    methods
        function this = RegWorker(caller_, volopt_)
            %REGISTER A constructor
            arguments
                caller_ (1,1)  Register
                volopt_ (1,12) table    % image volumes options
            end

            this.caller = caller_;
            this.volopt = volopt_;
            this.tform = cell(0,3);

            this.init_regopt(volopt_.cOrder, volopt_.slices);
        end

        function status = Run(this, sinfo_, regopt_, regfr_)
            % This function is the controller of registration
            arguments
                this
                sinfo_  (1,1)  struct    % source informaton struct
                regopt_ (1,1)  struct    % registration options
                regfr_  (1,:)  double {mustBePositive, mustBeInteger} 
            end

            % 1. initialize task manager
            this.taskmgr = TaskManager(this.volopt, regopt_)
        end

        function [status, path] = Save(this, fname, mdata, cidx)
            % This function save aligned movie
            arguments
                this
                fname   (1,1) string    % full file name
                mdata   (1,1) struct    % meta data struct
                cidx    (1,1) double {mustBeMember(cidx, [1,2])} = 1
            end

            data_export = squeeze(this.mov_aligned(:,:,cidx,:,:));

            % call the external function: savefile
            [status, path] = savefile(data_export,[],...
                'filename', fname, ...
                'metadata', mdata);

        end

        function status = Stop(this)
            this.state = matlab.lang.OnOffSwitchState("off");

            status = this.STATUS_SUCCESS;
        end

        function status = SetRegOpt(this, regopt)

        end
    end

    methods(Access = private)
        function init_regopt(this, c_order, slices)
            % init the registration options
            % init common parameters
            % ============ default registration type ==================
            this.regopt.RegType = "auto";

            % ============ automatic registration parameters ===========
            this.regopt.volfixed = struct("G", ["mean", "1"],...
                "L",  ["mean", "1"]);
            this.regopt.RegMode = "none";   % init as none, for no mode
            this.regopt.RegModal = "monomodal";
            
            this.regopt.AutoTform = "translation";
            this.regopt.InitStep = 0.1;
            this.regopt.MinStep = 1e-4;
            this.regopt.MaxIterN_Rigid = 50;
            this.regopt.MaxIterN_NonRigid = 100;
            this.regopt.IterCoeff = 0.5;
            this.regopt.AFS = 1.0;
            this.regopt.VPL_Rigid = floor(log(slices)/log(4))+1;
            this.regopt.VPL_NonRigid = this.regopt.VPL_Rigid;
            this.regopt.Interp_Rigid = "cubic";
            this.regopt.Interp_NonRigid = "cubic";

            % =========== manual registration parameters ==============
            this.regopt.ManualTform = "nonreflectivesimilarity";
            this.regopt.Degree = 3;
            this.regopt.Nlwm = 12;
            this.regopt.ManualInterp = "cubic";
            this.regopt.ManualEdgeSmooth = true;

            % init specific parameters
            switch numel(c_order)
                case 1
                    % ==== automatic registration parameters (1 channel) ====
                    this.regopt.RL = "ORN";
                    this.regopt.GA = 2.0;
                    this.regopt.InsTh = 2.0;
                    this.regopt.ScaleTh = 2.0;
                    this.regopt.GridUnit = "auto";
                case 2
                    % ==== automatic registration parameters (2 channels) ====
                    this.regopt.AutoContrast = false;
                    this.regopt.CompAcc = 1024;
                    this.regopt.CoRegC = 29;
                    this.regopt.DS = "auto";
                    this.regopt.Chla = "r";
                    this.regopt.Chls = "g";
                otherwise
            end

            % ============ long-term registration parameters =============
            this.regopt.KeyFrames = "auto";
            this.regopt.IntensityThreshold = 97;
            this.regopt.ScaleThreshold = 3;
            this.regopt.DS_PointCloud = "GA";
            this.regopt.DS_Voxel = "2X2";
            this.regopt.DS_PointCloud_Param = 2;
            this.regopt.OutlierRatio = 0.1;
            this.regopt.MaxIterN_PointCloud = 50;
            this.regopt.ErrorLimit = 1e-5;
            this.regopt.InitStep_LongTerm = 6.25e-2;
            this.regopt.MinStep_LongTerm = 1e-5;
            this.regopt.MaxIterN_Voxel = [100,50,25];
            this.regopt.IterCoeff_LongTerm = 0.5;
            this.regopt.AutoContrast_LongTerm = true;
            this.regopt.CompAcc_LongTerm = 1024;
            this.regopt.AFS_LongTerm = 1.0;
        end

        function status = align_one_vf(this)
            % This function align one virtual frame(time slice)

        end
    end
end

