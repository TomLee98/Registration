classdef TaskManager < handle
    %TASKMANAGER This class is task manager object defination
    % 
    
    properties(SetAccess=private, GetAccess=private)
        regopt      % 1-by-1 struct, the registration options
        volopt      % 1-by-12 table, the volume series options
        regfn       % 1-by-1 positive integer, registration frames number
        rcf         % 1-by-1 struct, the content of resource communication file
                    % note that RCF should be a public protocal, shared
                    % with different apps running
        taskqueue   % 1-by-1 mQueue, queue with Task object
    end

    properties(Hidden, Access=private)
        sfolder
        fname
        nmax

        caller
    end

    properties(Hidden, Constant)
        % STATUS DEFINATION
        STATUS_SUCCESS = 0
        STATUS_READ_FAIL = -1
        STATUS_WRITE_FAIL = -2
        STATUS_OTHER_PLATFORM = -3
        STATUS_NO_RESOURCE = -4
        STATUS_TOO_MANY_ALLOC = -5

        FIRST_ALLOC_COEFF = 0.5
        MAX_USER_ALLOC = 6

        % RCF STRUCT FIELD DEFINATION
        RCF_FIELD_NAMES = {'user_id', 'submit_time', 'time_use', ...
            'workers_num', 'is_alloc', 'is_gpu'}

        % GLOBAL REGISTRATION HIDDEN PARAMETERS
        IMREGTFORM_TRS_PPS = 2.60E5     % pixels per second
        IMWARP_TRS_CUBIC_PPS = 3.26E5
        IMWARP_TRS_LINEAR_CUBIC_RATIO = [3, 1];      
        CODENOISE_TRS_PPS = 1.56E6
        IMREGTFORM_TRS_RID_AFF_RATIO = [3,2,1]     % processing speed ratio: imregtform
        IMWARP_TRS_CUBIC_RID_AFF_RATIO  = [6,7,8]  % processing speed ratio: imwarp
        CODENOISE_TRS_RID_AFF_RATIO = [1,0,0]      % ~

        % LOCAL REGISTRATION HIDDEN PARAMETERS
        IMWARP_DEMON_TUSE = 10                      % typical using, actually can be omitted
        IMREGDEMONS_ACF_1_PPS = 1.04E5              % pixels per second
    end

    methods(Access = public)
        function this = TaskManager(volopt_, regopt_, regfn_, srvhost_)
            %TASKMANAGER A constructor
            arguments
                volopt_     (1,:)  table            % the movie information table
                regopt_     (1,1)  struct           % the registration options struct
                regfn_      (1,1)  double {mustBePositive, mustBeInteger}
                srvhost_    (1,1)  string = "silab" % server host name
            end

            this.volopt = volopt_;
            this.regopt = regopt_;
            this.regfn = regfn_;
            [~, this.nmax] = GetRunningMode();

            [~, host] = system("hostname");
            host = string(host).extractBefore(newline);
            if (this.nmax > 8) && (host == srvhost_)
                if isunix()
                    sfpath = '/data/.rcfs';
                    status = 0;             % avoid users permission
                elseif ispc()
                    sfpath = 'C:\User\Public\rcfs';
                    status = mkdir(sfpath);
                end

                if status == 0
                    throw(MException("TaskManager:mkdirFailed", ...
                        "Make dir failed."));
                end
                this.sfolder = string(sfpath);
                if ispc()
                    % hide the folder
                    fileattrib(this.sfolder, "+h");
                end
            else
                this.sfolder = [];  % no matter the sfolder
            end

            this.fname = [];
        end

        function SetCaller(this, caller_)
            arguments
                this
                caller_ (1,1) Register
            end
            this.caller = caller_;
        end
        
        function [status, n] = Allocate(this, wtype)
            % this function return the feasible number of parallel workers
            % and generate the allocate record
            % input:
            %   - wtype: 1-by-1 string, worker type, can be "cpu" or "gpu"
            % output:
            %   - status: indicate the function return state
            %   - n: 1-by-x, the suggestion maximum number of workers, x
            %   can be 1 or 2, 1 for 'global-only', 'local-only','longterm'
            %   and 2 for 'mix'

            arguments
                this
                wtype (1,1) string {mustBeMember(wtype, ["cpu", "gpu"])}...
                     = "cpu"
            end

            if this.nmax <= 8
                n = this.nmax;
                status = this.STATUS_OTHER_PLATFORM;
                return;
            end

            % On Unix computational platfrom
            % 0. get the spacial limitation of cpu/gpu number
            switch wtype
                case "cpu"
                    is_bigfile = (this.regopt.BigFile == "on");
                    wnmax_sp = GetCpuWorkersMaxN(this.volopt, ...
                                                 is_bigfile);
                case "gpu"
                    wnmax_sp =GetGPUWorkersMaxN(this.volopt);
                otherwise
            end

            % 1. use rlmap to calculate the computatioanl resource needed
            cr = this.rlmap();

            % hold until there is resource for calculating
            while true
                % 2. read the RCF files to get allocated information
                [status, info] = this.read();

                if status == this.STATUS_SUCCESS
                    % 3. compare info and generate rcf struct
                    status = genrcf(this, info, cr, wnmax_sp, wtype);

                    if status == this.STATUS_SUCCESS
                        break;
                    end
                end
                fprintf("TaskManager:waiting for allocating computation resource...\n");

                pause(3);
            end

            % 4. write RCF to folder, return n
            status = this.write();
            n = this.rcf.workers_num;

            fprintf("TaskManager:Allocating computation resource success.\n");
        end

        function status = Free(this)
            % this function remove the allocate record, which must be called
            % after registration for releasing resource

            if ~isempty(this.fname) && isfile(this.fname)
                delete(this.fname);
                this.fname = [];
            end

            status = this.STATUS_SUCCESS;
        end
    end

    methods(Access = private)
        function cr = rlmap(this)
            % this function generate the resource level mapping:
            % (volopt, regopt) -> positive real number (normalized calculation time)

            % unit: *100 seconds
            cr = (this.codenoise_time_use() ...
                + this.imregtform_time_use() ...
                + this.imwarp_time_use() ...
                + this.imregdemons_time_use()) / 100;
        end

        function t_use = imregtform_time_use(this)
            switch this.regopt.RegType
                case "auto"
                    if ismember(this.regopt.RegMode, ["global-only","mix"])
                        [~, rp] = ismember(this.regopt.AutoTform, ...
                            ["translation","rigid","similarity","affine"]);
                        if rp == 4, rp =3; end
                        procsp = this.IMREGTFORM_TRS_PPS ...
                            *this.IMREGTFORM_TRS_RID_AFF_RATIO(rp);
                        t_use = this.pixels_calc_num() / procsp;
                    else
                        t_use = 0;
                    end
                case "longterm"
            end
        end

        function t_use = codenoise_time_use(this)
            if ismember(this.regopt.RegMode, ["global-only","mix"]) ...
                    && (this.regopt.AutoTform == "translation")
                t_use = this.pixels_calc_num() ...
                    / this.CODENOISE_TRS_PPS;
            else
                t_use = 0;
            end
        end

        function t_use = imwarp_time_use(this)
            [~, rp_type] = ismember(this.regopt.AutoTform, ...
                ["translation","rigid","similarity","affine"]);
            if rp_type == 4, rp_type =3; end

            switch this.regopt.RegMode
                case "global-only"
                    [~, rp_alg] = ismember(this.regopt.Interp_Rigid, ["linear","cubic"]);
                    procsp = this.IMWARP_TRS_CUBIC_PPS ...
                        *this.IMWARP_TRS_LINEAR_CUBIC_RATIO(rp_alg) ...
                        *this.IMREGTFORM_TRS_RID_AFF_RATIO(rp_type);
                    t_use = this.pixels_calc_num() / procsp;
                case "local-only"
                    t_use = this.IMWARP_DEMON_TUSE;
                case "mix"
                    [~, rp_alg] = ismember(this.regopt.Interp_Rigid, ["linear","cubic"]);
                    procsp = this.IMWARP_TRS_CUBIC_PPS ...
                        *this.IMWARP_TRS_LINEAR_CUBIC_RATIO(rp_alg) ...
                        *this.IMREGTFORM_TRS_RID_AFF_RATIO(rp_type);
                    t_use = this.pixels_calc_num() / procsp ...
                        + this.IMWARP_DEMON_TUSE;
            end
        end

        function t_use = imregdemons_time_use(this)
            switch this.regopt.RegMode
                case {'local-only', 'mix'}
                    pn =  this.volopt.width*this.volopt.height*this.volopt.slices*this.regfn;
                    t_use = pn / this.IMREGDEMONS_ACF_1_PPS;
                    t_use = 0.5*t_use *(this.regopt.AFS+1);
                otherwise
                    t_use = 0;
            end
        end

        function pn = pixels_calc_num(this)
            if this.volopt.channels == 1

            elseif this.volopt.channels == 2
                if this.regopt.DS == "auto"
                    scale = 256./max([this.volopt.width, ...
                        this.volopt.height]);
                else
                    scale = str2double(string(this.regopt.DS).extractBefore(2));
                end

                pn =  this.volopt.width*this.volopt.height/scale.^2 ...
                    *this.volopt.slices ...
                    *this.regfn;
            end
            
        end

        function status = write(this)
            % Write This function write resource communication file (RCF)
            % to sfolder, name as userid_datetime.xml

            fname_ = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];
            fname_ = [char(fname_(randperm(18))), '.xml'];
            this.fname = string([this.sfolder.char(), filesep, fname_]);

            try
                writestruct(this.rcf, this.fname, "FileType","xml");
            catch ME
                throwAsCaller(ME);
            end

            status = this.STATUS_SUCCESS;
        end

        function [status, info] = read(this)
            % Read This function read the total resource communicatio
            % files (RCF) from sfolder
            % output:
            %   - status: indicate the function return state
            %   - info: struct, with total acruired resource and total cr

            info = [];

            % read all *.xml files in sfolder
            xml_files = struct2table(dir(this.sfolder));
            [~,~,ext] = fileparts(xml_files.name);
            ext = string(ext);
            xml_files = xml_files.name(ext.contains(".xml"));

            for n = 1:numel(xml_files)
                % struct with fields:
                % 'submit_time', 'task_use', 'workers_num', 'is_alloc'
                s = readstruct([this.sfolder.char(), filesep, xml_files{n}], ...
                    "FileType","xml");
                validatercf(s);

                info = [info; s]; %#ok<AGROW>
            end

            status = this.STATUS_SUCCESS;

            function validatercf(s)
                % s should be a struct
                if ~isempty(setxor(fieldnames(s), this.RCF_FIELD_NAMES))
                    throw(MException("TaskManager:read:invalidRCFFile", ...
                        "Some *.xml files may be damaged."));
                end

                validateattributes(s.user_id, "string", "scalar");
                validateattributes(s.submit_time, "datetime","scalar");
                validateattributes(s.time_use, "double", "scalar");
                validateattributes(s.workers_num, "double", {'scalar', 'integer','positive'});
                validateattributes(s.is_alloc, "logical", "scalar");
                validateattributes(s.is_gpu, "logical", "scalar");
            end
        end

        function status = genrcf(this, info, cr, n_ref, wtype)
            % This function generate rcf struct
            % input:
            %   - info: n-by-1 struct, others rcf
            %   - cr: 1-by-1 double, this process cost
            %   - n_ref: 1-by-1 positive integer, the advised maximum
            %            allocating workers number

            if ispc()
                [~, uid] = system("whoami /logonid");
            elseif isunix()
                [~, uid] = system("id -u");
            end
            uid = string(uid).extractBefore(newline);
            this.rcf.user_id = uid;

            workers_alloc = 0;
            cr_alloc = 0;
            user_alloc_num = 0;
            for k = 1:numel(info)
                if info(k).is_alloc == true
                    workers_alloc = workers_alloc + info(k).workers_num;
                    cr_alloc = cr_alloc + info(k).time_use;
                    user_alloc_num = user_alloc_num + 1;
                end
            end

            if user_alloc_num == this.MAX_USER_ALLOC
                status = this.STATUS_TOO_MANY_ALLOC;
                return;
            end

            if workers_alloc == 0
                % first allocate resource
                workers_to_alloc = round(this.nmax * this.FIRST_ALLOC_COEFF);
            else
                % compare the time_use and allocated workers, prop.
                workers_to_alloc = round(cr / cr_alloc*workers_alloc);

                % allocate the left workers
                workers_to_alloc = min(workers_to_alloc, ...
                    this.nmax - workers_alloc);
            end
            this.rcf.workers_num = min([workers_to_alloc, n_ref, this.volopt.frames]);

            this.rcf.submit_time = datetime();

            this.rcf.time_use = cr;

            this.rcf.is_gpu = (wtype == "gpu");     % ? what effect?

            if this.rcf.workers_num > 0
                this.rcf.is_alloc = true;
                status = this.STATUS_SUCCESS;
            else
                this.rcf.is_alloc = false;
                status = this.STATUS_NO_RESOURCE;
            end
        end
    end
end

