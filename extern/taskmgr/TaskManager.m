classdef TaskManager < handle
    %TASKMANAGER This class is task manager object defination
    
    properties(SetAccess=private, GetAccess=private)
        regopt      % 1-by-1 struct, the registration options
        volopt      % 1-by-12 table, the volume series options
        rcf         % 1-by-1 struct, the content of resource communication file
    end

    properties(Hidden, Access=private)
        sfolder
        nmax

        caller
    end

    properties(Hidden, Constant)
        STATUS_SUCCESS = 0
        STATUS_READ_FAIL = -1
        STATUS_WRITE_FAIL = -2
        STATUS_OTHER_PLATFORM = -3

        % GLOBAL REGISTRATION HIDDEN PARAMETERS
        IMREGTFORM_TRS_PPS = 2.60E5     % pixels per second
        IMWARP_TRS_PPS = 3.26E5         
        CODENOISE_TRS_PPS = 1.56E6
        IMREGTFORM_TRS_RID_AFF_RATIO = [3,2,1]     % processing speed ratio: imregtform
        IMWARP_TRS_RID_AFF_RATIO  = [6,7,8]        % processing speed ratio: imwarp
        CODENOISE_TRS_RID_AFF_RATIO = [1,0,0]      % ~

        % LOCAL REGISTRATION HIDDEN PARAMETERS
        
    end
    
    methods(Access = public)
        function this = TaskManager(sfolder_, regopt_, volopt_, nmax_, hostname_)
            %TASKMANAGER A constructor
            arguments
                sfolder_ (1,1)  string    % the RCF folder
                regopt_  (1,1)  struct    % the registration options struct
                volopt_  (1,12) table     % the movie information table
                nmax_    (1,1)  double {mustBePositive, mustBeInteger} = 48
                hostname_(1,1)  string = "silab"
            end

            this.volopt = volopt_;
            this.regopt = regopt_;
            this.nmax = nmax_;

            [~, host] = system("hostname");
            host = string(host).extractBefore(newline);
            if (nmax_ > 8) && (host == hostname_)
                if isfolder(sfolder_)
                    this.sfolder = sfolder_;
                else
                    warning("TaskManager:invalidFolder", "The folder is invalid, " + ...
                        "shared folder will be generated automatically.");
                    if isunix()
                        sfpath = '/data/rcfs';
                    elseif ispc()
                        sfpath = 'C:\User\Public\rcfs';
                    end
                    status = mkdir(sfpath);
                    if status == 0
                        throw(MException("TaskManager:mkdirFailed", ...
                            "Make dir failed."));
                    end
                    this.sfolder = string(sfpath);
                end
            else
                this.sfolder = "";  % no matter the sfolder
            end
        end

        function SetCaller(this, caller_)
            arguments
                this
                caller_ (1,1) Reg3D_GUI
            end
            this.caller = caller_;
        end
        
        function [status, n] = Allocate(this)
            % this function return the feasible number of parallel workers
            % output:
            %   - status: indicate the function return state
            %   - n: the suggestion maximum number of workers

            if this.nmax <= 8
                n = this.nmax;
                status = this.STATUS_OTHER_PLATFORM;
                return;
            end

            % On Unix computational platfrom

            % 1. use rlmap to calculate the computatioanl resource needed
            cr = this.rlmap();

            % 2. read the RCF files to get allocated information
            [status, info] = this.read();

            if status == this.STATUS_SUCCESS
                % 3. compare info and generate rcf struct
                genRCF(this, info, cr);

                % 4. write RCF to folder, return n
                status = this.write();

                n = this.rcf.worker_n;
            end
        end
    end

    methods(Access = private)
        function cr = rlmap(this)
            % this function generate the resource level mapping:
            % (volopt, regopt) -> positive real number

            % 1. registration speed relationship?
            % (1) normal registration mode
            % [1] change frames number: 
            % [2] change volume size
            % [3] 

        end

        function status = write(this)
            % Write This function write resource communication file (RCF)
            % to sfolder, name as userid_datetime.xml

            this.rcf.is_alloc = true;
        end

        function [status, info] = read(this)
            % Read This function read the total resource communicatio
            % files (RCF) from sfolder
            % output:
            %   - status: indicate the function return state
            %   - info: struct, with total acruired resource and total cr

        end

        function status = acquire(this)
            %Acquire This function edit RCF and save the require recording

        end

        function status = genRCF(this, info, cr)
            this.rcf
        end
    end
end

