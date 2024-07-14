classdef regchain < handle
    %REGCHAIN This class defines registration chain object, and we can do
    %some operation on this object
    % generate the registration chain

    properties(Access=private, Hidden, Constant)
        REGRESSION_SIGMA_N = 2
    end

    properties(GetAccess=public, Dependent)
        isValid         % logical
        Template        % regtmpl object
        ChainInfo       % struct, with {data, info}
    end
    
    properties(Access=private, Hidden)
        movsrc      % 1-by-1 regmov, the movie source
        tbl_chain   % k-by-4 table, with {kfidx, fscope, tfmat, kfvol}
        regtpl      % 1-by-5 cell, {absolute_index, channel_index, 3D array, mass_center, mass_shift}
        opt_chain   % 1-by-1 struct, with registration chain options
        empty_flag  % 1-by-1 logical, indicate the empty status
        valid_flag  % 1-by-1 logical, indicate the validation status
    end

    methods
        function r = get.isValid(this)
            r = this.valid_flag;
        end

        function r = get.Template(this)
            r = regtmpl(this.movsrc, struct("Sampling",["mean", string(num2str(this.regtpl{1}))], ...
                                            "Channel", this.opt_chain.sc));
        end

        function r = get.ChainInfo(this)
            r = struct("data", this.tbl_chain, ...
                       "info", {this.regtpl});
        end
    end

    methods(Access=public)
        function this = regchain(movsrc_, movtmpl_, opts_, isemp_)
            %REGCHAIN A constructor
            arguments
                movsrc_     (1,1)   regmov
                movtmpl_    (1,1)   regtmpl
                opts_       (1,1)   struct
                isemp_      (1,1)   logical = false
            end

            this.movsrc = movsrc_;

            % get volume and index information
            fidef = movtmpl_.FixDef;
            this.regtpl = {median(str2num(fidef.Sampling(2))), ...
                           find(fidef.Channel==movsrc_.MetaData.cOrder), ...
                           movtmpl_.RefVol, ...
                           [0,0,0,0], ...
                           nan(movsrc_.MetaData.frames,1)}; %#ok<ST2NM>

            % reorganize the chain options
            this.opt_chain = struct("template_auto", opts_.AutoTemplate, ...
                                    "keyframe_auto", opts_.AutoKeyframe, ...
                                    "keyframes",     opts_.Keyframes, ...
                                    "tgrid_minmax",  opts_.TGridMinMax, ...
                                    "sc",            opts_.SC, ...
                                    "fc",            opts_.FC, ...
                                    "iter_coeff",    opts_.IterCoeff, ...
                                    "itn_max",       opts_.MaxIterN, ...
                                    "step_max",      opts_.MaxStep, ...
                                    "step_min",      opts_.MinStep, ...
                                    "zopt_shift_max",opts_.MaxZOptShift, ...
                                    "dfsize",        opts_.DilateFilter);

            this.empty_flag = isemp_;
            this.valid_flag = false;
        end

        function generate(this)
            if isempty(this)
                throw(MException("regchain:emptyObject", "Generation is invalid " + ...
                    "on empty object."))
            end

            % This is core function, which will generate a registration
            % chain and store in tbl_chain

            %   1 select the template
            %       1.1   given by user
            %       1.2   autometic (nearest the 'mass center' and 'time center' in S-T)
            if this.opt_chain.template_auto == true     % overwrite the template
                % get all volumes mass center
                mc = this.movsrc.MC(:,:,this.regtpl{2});    % t-by-3

                % nearest the spatial-temporal median frame
                rms = sqrt(mean(var(mc(:,1:3), 0, 1)));
                t_cpb = (1:size(mc,1))';
                t_cpb = t_cpb*rms/sqrt(numel(t_cpb)*(numel(t_cpb)+1)/12);
                t_cpb = t_cpb - mean(t_cpb);    % zero mean

                % combine spatial and temporal information
                mc_cpb = [mc-mean(mc,1), t_cpb];
                mc = [mc, (1:size(mc,1))'];

                [~, pmc_st] = min(std(mc_cpb, 0, 2));

                % update regtpl
                this.regtpl{1} = pmc_st;
                this.regtpl{3} = grv(this.movsrc, ["none", string(pmc_st)], ...
                    this.opt_chain.sc);
                this.regtpl{4} = mc(pmc_st, :);
                this.regtpl{5} = sqrt(sum(mc_cpb(:,1:3).^2, 2));
            end

            %   2 select the keyframes
            %       2.1   given by user
            %       2.2   autometic
            %           2.2(a)  generate partitions s.t. displacement variations
            %                     in partitions are similar, and the mean of variations
            %                     is lower than estimation with given registration
            %                     parameters (iteration, partitions number range)
            %           2.2(b)  select the keyframe in each partition, which obey the
            %                     same rule for template generation(self-similarity)
            if this.opt_chain.keyframe_auto == true
                if ~exist("mc", "var")
                    % generate the mass center
                    mc = this.movsrc.MC;    % t-by-3-by-c
                    mc = [mc(:,:,this.regtpl{2}), (1:size(mc,1))'];
                end

                % compress the motion as results estimation after coarse
                % registration
                NormCoeff = [this.movsrc.MetaData.width, this.movsrc.MetaData.height, ...
                    this.movsrc.MetaData.slices, this.movsrc.MetaData.frames];
                mc_f = (mc - this.regtpl{4})./NormCoeff;     % normlized to [-1, 1]
                mc_fd = sqrt(sum(mc_f(:,1:3).^2, 2));

                % auto refine the partition
                % generate the initial partition:
                NP = this.opt_chain.tgrid_minmax(1);
                NPm = this.opt_chain.tgrid_minmax(2);
                fparts = ceil(linspace(1, this.movsrc.MetaData.frames+1, NP)); % LCRO
                % assume that MC variation could be explained by
                % translation and shear with same weight
                % then var(total) = 2 * var(translation)
                dth = sqrt(this.opt_chain.step_max^2/(1-this.opt_chain.iter_coeff^2)*2);
                pdv = partial_spd(fparts);  % partial speed in 1 time step, <=> unit displacement

                % decrease the variation until it's satisfied affine 
                % registration or come to partition number upper bound
                while any(pdv > dth) && (NP < NPm)
                    % insert a new partition into old partition with the
                    % maximum variation
                    [~, vloc] = max(pdv);

                    % find the best seperation point, fine sampling at the
                    % high speed point
                    mc_sfd = mc_fd(fparts(vloc):fparts(vloc+1)-1);
                    [~, dvloc] = max(abs(diff(mc_sfd)));    % 1-order diff, left shift 1 bit

                    % update fparts: insert
                    fparts = [fparts(1:vloc), fparts(vloc)+dvloc, fparts(vloc+1:end)];

                    NP = NP + 1;
                    pdv = partial_spd(fparts);
                end

                % select keyframes in partition: fparts
                kfidx = nan(numel(fparts)-1, 1);
                fscope = nan(numel(fparts)-1, 2);
                for n = 1:numel(kfidx)
                    mc_sfd = mc_fd(fparts(n):fparts(n+1)-1);
                    [~, dvloc] = min(abs(mc_sfd-median(mc_sfd)));
                    kfidx(n) = fparts(n) + dvloc - 1;
                    fscope(n, :) = [fparts(n), fparts(n+1)-1];
                end
            else
                kfidx = this.opt_chain.keyframes;
                
                % find the symmetric fscope to contain kfidx
                fscope = nan(numel(kfidx), 2);
                pfscope = nan(1, numel(kfidx)-1);
                for n = 1:numel(kfidx)-1
                    pfscope(n) = round((kfidx(n)+kfidx(n+1))/2);
                end
                pfscope = [1, pfscope, this.movsrc.MetaData.frames+1];
                for n = 1:numel(kfidx)
                    fscope(n, :) = [pfscope(n), pfscope(n+1)-1];
                end
            end

            %   3 use phase correlation registration for rigid body transformation as
            %       initial value, then do intensity based affine transformation
            %       * use two process track previous and posterior temporal parts
            [tfmat, kfvol] = align_chain_node(this, kfidx);

            % generate the chain table
            this.tbl_chain = table(kfidx, fscope, tfmat, kfvol);

            this.valid_flag = true;

            function pv = partial_spd(part)
                pv = nan(1, numel(part)-1); % 1,2,...,n-1
                for k = 1:numel(pv)
                    % cover 95% displacement field unit under exponential
                    % distribution
                    pv(k) = 3*std(diff(mc_fd(part(k):part(k+1)-1)));
                end
            end

        end

        function [tf, kv] = acquire(this, fidx)
            arguments
                this (1,1)
                fidx (1,1)  double  {mustBeInteger, mustBePositive}
            end

            if ~isempty(this) && this.isValid
                % This function return the best keyframe and pre-aligned
                % transformation matrix
                loc = (fidx>=this.tbl_chain.fscope(:,1) ...
                    & fidx<=this.tbl_chain.fscope(:,2));
                if any(loc)     % infect, trange has only one partition solution
                    tf = this.tbl_chain.tfmat(loc);
                    kv = this.opt_chain.kfvol(loc);
                else
                    % could never come here
                    tf = eye(size(this.tbl_chain.kfvol(1))); % identity matrix
                    kv = this.regtpl{3};                     % template volume
                end
            else
                throw(MException("regchain:emptyObject", "Requirement is invalid " + ...
                    "on empty object."));
            end
        end

        function [tfs, kvs, ids] = vacquire(this, vfidx)
            arguments
                this    (1,1)
                vfidx   double  {mustBeVector, mustBeInteger, mustBePositive}
            end

            if ~isempty(this) && this.isValid
                tfs = cell(size(this.tbl_chain, 1), 1);
                kvs = cell(size(tfs));
                ids = nan(size(vfidx));
                for k = 1:numel(vfidx)
                    loc = (vfidx(k)>=this.tbl_chain.fscope(:,1) ...
                        & vfidx(k)<=this.tbl_chain.fscope(:,2));
                    ids(k) = find(loc);
                    tfs{ids(k)} = this.tbl_chain.tfmat{loc};
                    kvs{ids(k)} = this.tbl_chain.kfvol{loc};
                end
            else
                 throw(MException("regchain:emptyObject", "Requirement is invalid " + ...
                    "on empty object."));
            end
        end

        function tf = isempty(this)
            tf = this.empty_flag;
        end
    end

    methods(Access = private)
        function [tf, kv] = align_chain_node(this, kfidx)
            % This function use 'phase correlation' registration as coarse
            % registration and affine registration as fine tuning
            tf = cell(numel(kfidx), 1);
            kv = cell(numel(kfidx), 1);

            fmode_pre_inv = ["none", string(flipud(kfidx(kfidx<=this.regtpl{1}))).join(",")];   % reverse time order
            fmode_pre = ["none", string(kfidx(kfidx<=this.regtpl{1})).join(",")];   % reverse time order
            n_pre = sum(kfidx<=this.regtpl{1});
            fmode_post = ["none", string(kfidx(kfidx>this.regtpl{1})).join(",")];
            n_post = sum(kfidx>this.regtpl{1});
            
            tfc = {repmat(eye(4),1,1,n_pre+1), repmat(eye(4),1,1,n_post+1)};    % pre, post

            avol_sc_pre_inv = grv(this.movsrc, fmode_pre_inv, this.opt_chain.sc);
            avol_sc_pre_inv = cat(4, this.regtpl{3}, avol_sc_pre_inv);
            avol_sc_pre = grv(this.movsrc, fmode_pre, this.opt_chain.sc);
            avol_sc_post = grv(this.movsrc, fmode_post, this.opt_chain.sc);
            avol_sc_post = cat(4, this.regtpl{3}, avol_sc_post);

            avol = {avol_sc_pre_inv, avol_sc_post};

            % unzip the regstration options
            iter_coeff = this.opt_chain.iter_coeff;
            itn_max = this.opt_chain.itn_max;
            step_max = this.opt_chain.step_max;
            step_min = this.opt_chain.step_min;
            zs_max = this.opt_chain.zopt_shift_max;
            dfsize = this.opt_chain.dfsize;
            rs = [this.movsrc.MetaData.xRes, this.movsrc.MetaData.yRes, ...
                this.movsrc.MetaData.zRes];
            rref = imref3d(size(avol_sc_pre_inv, 1:3), rs(1),rs(2),rs(3)); % assume vol: x,y,z,t
            vpl = fix(log(size(avol_sc_pre_inv,3))/log(4))+1;

            % alignment slide the chain
            % using parpool for maximum 2 fold speed up
            parobj = parpool("Processes",2, "SpmdEnabled",false);

            pbar = parfor_wait(numel(kfidx), 'Waitbar', true, ...
                'Name','Reg-Chain Generating', 'ReportInterval', 1);

            parfor worker = 1:2     % two workers work on different time lines
                N = {n_pre, n_post};

                % pre-loading data
                vol_fix = avol{worker}(:,:,:,1);
                
                % set up the registration optimizer options
                [optimizer, metric] = imregconfig("monomodal");
                optimizer.MaximumStepLength = step_max;
                optimizer.MinimumStepLength = step_min;
                optimizer.RelaxationFactor = iter_coeff;
                optimizer.MaximumIterations = itn_max;

                for nw = 2:N{worker}+1
                    fvol = avol{worker}(:,:,:,nw-1);
                    mvol = avol{worker}(:,:,:,nw);

                    % preprocessing for robust registration
                    fvol_proc = preproc_tc(fvol, dfsize);
                    mvol_proc = preproc_tc(mvol, dfsize);

                    [ptf, ~] = imregcoarse(mvol_proc, fvol_proc, rs, zs_max);

                    % do fine registration: imregtform base on intensity
                    ptf = imregtform(mvol, rref, fvol, rref, "affine", ...
                        optimizer, metric, "PyramidLevels",vpl, "InitialTransformation",ptf);

                    % multiply with pre-aligned transformation as
                    % initial matrix
                    % left multiply is wrong when MATLAB version < R2022b
                    ptf.A = tfc{worker}(:,:,nw-1)*ptf.A; %*v

                    % fine tuning: align to template frame
                    ptf = imregtform(mvol, rref, vol_fix, rref, "affine", ...
                        optimizer, metric, "PyramidLevels",vpl, "InitialTransformation",ptf);

                    tfc{worker}(:,:,nw) = ptf.A;

                    pbar.Send; %#ok<PFBNS>
                    pause(0.002);
                end
            end

            pbar.Destroy;

            % remove the header transformation
            tfc{1}(:,:,1) = [];     tfc{2}(:,:,1)=[];
            avol_sc_post(:,:,:,1) = [];
            
            % modify the output data
            for m = 1:n_pre
                tf{m} = tfc{1}(:,:,n_pre-m+1);
                kv{m} = avol_sc_pre(:,:,:,m);
            end
            for m = n_pre+1:numel(kfidx)
                tf{m} = tfc{2}(:,:,m-n_pre);
                kv{m} = avol_sc_post(:,:,:,m-n_pre);
            end

            delete(parobj);
        end
    end

    methods(Static, Hidden)
        function c = empty()
            % this function generates an empty registration chain
            opts = struct("AutoTemplate", true, ...
                          "AutoKeyframe", true, ...
                          "Keyframes",    1, ...
                          "TGridMinMax",  [2, 100], ...
                          "SC",           "r", ...
                          "FC",           "g", ...
                          "IterCoeff",    0.5, ...
                          "MaxIterN",     50, ...
                          "MaxStep",      0.0625, ...
                          "MinStep",      1e-7, ...
                          "MaxZOptShift", 2, ...
                          "DilateFilter", [3,115,1000]);

            c = regchain(regmov.empty(), ...
                         regtmpl.empty(), ...
                         opts, ...
                         true);
        end
    end
end

