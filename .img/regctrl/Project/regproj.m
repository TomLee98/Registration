classdef regproj < handle
    %REGPROJ This Class define Project object, which save data related to
    %a registration project
    
    properties (GetAccess = public, Dependent)
        LastSavedTime   % get/___,     stored, 1-by-1 datetime indicate the last saved time
        Seed            % get/___, not stored, 1-by-1 struct could be save as project
        IsReady         % get/___,     stored, 1-by-1 logical indicate if all data is ready
        IsSaved         % get/set,     stored, 1-by-1 logical indicate if project is saved
    end

    % project immutable variables:
    % set only once and will not change in lifetime
    properties(Access = {?regpm, ?Register})
        pcolor          (1,1)   struct  = struct("normal",zeros(3), ...
                                                 "fused", zeros(3))                 % pseudo RGB color definition, {normal, fused}
        proj_fname      (1,1)   string  = ""                                        % current project filename
    end

    % project active variables:
    % set according to user operations
    properties(Access = {?regpm, ?Register})
        %% variables could restore after loading
        cidx_cur        (1,:)   logical  {mustBeNumericOrLogical} = []              % 1-by-c logical array, the current channel index
        c_order         (1,:)   string  {mustBeMember(c_order, ["r","g","b"])} = [] % 1-by-c string, the channels order
        fixdef          (1,1)   struct  = struct("Sampling",["mean","1"], ...
                                                 "Channel","r");                    % the template volume description
        fidx_cur        (1,1)   double  {mustBeInteger, mustBeNonnegative} = 0      % the current frame index
        fps             (1,1)   struct  = struct("fold",1, "base",[])               % the play speed setting, {base, fold} 
        frames_reg      (1,2)   string  = ["1:end", "1"]                            % the registration frames, [auto, manual]
        img_cur         (1,1)   regmov  = regmov.empty()                            % current active image pointer
        img_fname       (1,1)   string  = ""                                        % record the processing raw image file name
        loaded_flag     (1,1)   logical = false                                     % flag for image data status, true for loaded
        opt_crop        (1,1)   struct  = struct("xy",[],"z",[],"f",[])             % the current ROI for image crop, field with xy, z and f(frame)
        opt_movout      (1,1)   struct  = struct()                                  % the movie exporting configuration
        opt_reg         (1,1)   regopt  = regopt("TCREG", "global")                 % registration options
        opt_seg         (1,1)   segopt  = segopt()                                  % segmentation options
        opt_sig         (1,1)   sigopt  = sigopt("MovingQuantile", "Constant")      % signal extraction options
        refvol          (1,1)   regtmpl = regtmpl.empty()                           % current registration template object
        saved_flag      (1,1)   logical = false                                     % 1-by-1 logical, flag for project saved status
        saved_time      (1,1)   datetime= datetime("tomorrow")                      % 1-by-1 datetime record project saved time
        sidx_cur        (1,1)   double  {mustBeInteger, mustBeNonnegative} = 0      % the current slice index
        
        %%  variables could restore after development
        dev_done        (1,1)       logical = false                                 % flag for if data is developed
        img_zproj       (1,:)       struct  = struct("cur",[], "ref", []);          % the z-projected images, need memory to restore
        mimg_cur        (:,:,:)     uint16  = []                                    % the coarse registration inner volume
    end
    
    methods
        %% Constructor
        function this = regproj()
            %REGPROJ A Constructor
        end
        
        %% Getter / Setter

        function r = get.LastSavedTime(this)
            if this.saved_flag == false && this.saved_time<= datetime("now")
                % reset
                this.saved_time = datetime("tomorrow");
            end

            % save not now
            r = this.saved_time;
        end

        function r = get.Seed(this)
            r = struct("cidx_cur",this.cidx_cur, "c_order",this.c_order, ...
                "fixdef",this.fixdef, "fidx_cur",this.fidx_cur, ...
                "fps",this.fps, "frames_reg",this.frames_reg, ...
                "img_cur",this.img_cur, "img_fname",this.img_fname, ...
                "loaded_flag",this.loaded_flag, "opt_crop",this.opt_crop, ...
                "opt_movout",this.opt_movout, "opt_reg",this.opt_reg, ...
                "opt_seg",this.opt_seg, "opt_sig",this.opt_sig, ...
                "pcolor",this.pcolor, "proj_fname",this.proj_fname, ...
                "refvol",this.refvol, "saved_flag",this.saved_flag, ...
                "saved_time",this.saved_time, "sidx_cur",this.sidx_cur);
        end

        function r = get.IsReady(this)
            r = this.dev_done;
        end

        function r = get.IsSaved(this)
            r = this.saved_flag;
        end

        function set.IsSaved(this, r)
            arguments
                this
                r       (1,1)   logical
            end

            this.saved_flag = r;
            if this.saved_flag == true
                this.saved_time = datetime("now");
            end
        end
    end

    methods (Access = public)
        function Restore(this, seed)
            arguments
                this
                seed    (1,1)   struct
            end

            % distribute stored varibles
            this.redist(seed);

            % develop other variables
            this.develop();
        end

        function Reset(this)
            pcolor_ = this.pcolor;
            proj_fname_ = this.proj_fname;

            % reset current handle
            this = regproj.empty();

            this.pcolor = pcolor_;
            this.proj_fname = proj_fname_;
        end

        function r = isempty(this)
            r = ~this.loaded_flag;
        end

        function delete(this)
            % free link
            delete(this.img_cur);
            delete(this.refvol);

            % ~
        end
    end

    methods (Access = private)
        function redist(this, seed)
            this.cidx_cur = seed.cidx_cur;
            this.c_order = seed.c_order;
            this.fixdef = seed.fixdef;
            this.fidx_cur = seed.fidx_cur;
            this.fps = seed.fps;
            this.frames_reg = seed.frames_reg;
            this.img_cur = seed.img_cur;
            this.img_fname = seed.img_fname;
            this.loaded_flag = seed.loaded_flag;
            this.opt_crop = seed.opt_crop;
            this.opt_movout = seed.opt_movout;
            this.opt_reg = seed.opt_reg;
            this.opt_seg = seed.opt_seg;
            this.opt_sig = seed.opt_sig;
            this.pcolor = seed.pcolor;
            this.proj_fname = seed.proj_fname;
            this.refvol = seed.refvol;
            this.saved_flag = seed.saved_flag;
            this.saved_time = seed.saved_time;
            this.sidx_cur = seed.sidx_cur;
        end

        function develop(this)
            % develop img_zproj variable
            if this.img_cur.IsProjected
                this.img_zproj.cur = this.img_cur.MovieZProj;
                this.img_zproj.ref = Projection(this.refvol.RefVol, ...
                    this.img_cur.ZProjMethod, 3);
            end

            % develop mimg_cur
            vol = this.img_cur.Movie(:, :, this.cidx_cur, :, this.fidx_cur);
            vol = squeeze(vol);
            this.mimg_cur = GenPreProcVol(vol, this.opt_reg);

            this.dev_done = true;
        end
    end

    methods(Static)
        function r = empty()
            r = regproj();
        end
    end
end

