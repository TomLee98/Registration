classdef Estimator < handle
    %SIGNALEXTRACTOR This class define a signal extractor object

    properties(Constant, Hidden)

    end

    properties(GetAccess=public, Dependent)
        Activities      % get,  variable,   s-by-t double array
        Baseline        % get,  variable,   s-by-t double array
        Noise           % get,  variable,   s-by-t double array
        CaImAnData      % get,  variable,   1-by-1 struct, with field {A, C, b, f}
        Significant     % get,  variable,   s-by-1 logical array
    end
    
    properties(Access=private, Hidden)
        fl              % s-by-t double array, row for components, colume for time step
        dff             % s-by-t double array, row for components, colume for time step
        bl              % s-by-t double array, row for components, colume for time step
        noise           % s-by-t double array, row for components, colume for time step
        cidx            % 1-by-1 positive integer, indicate the channel index
        opts            % 1-by-1 sigopt object

        A_caiman        % m*n*p-by-s sparse double label matrix, s for number of components
        C_caiman        % s-by-t double matrix, indicates weighted components activities
        b_caiman        % m*n*p-by-r sparse array, r for rank of background
        f_caiman        % r-by-t double matrix, indicates weighted background activities
    end

    properties(SetAccess=immutable, Hidden)
        image_src       % 1-by-1 regmov object, the calcium image data 
        mask            % 1-by-1 regmask object
        caller          % 1-by-1 Register object
    end
    
    methods
        function this = Estimator(src_, mask_, caller_)
            arguments
                src_    (1,1)   regmov
                mask_   (1,1)   regmask
                caller_ (1,1)   Register
            end

            this.image_src = src_;
            this.mask = mask_;
            this.fl = nan(this.mask.NumComps, this.image_src.MetaData.frames);

            this.caller = caller_;
        end
        
        function fit(this, opts_, fc_, comps_, withcm_)
            % This function fit the signal model
            % DO NOT USE auto_parpool in this block
            arguments
                this
                opts_   (1,1)   sigopt
                fc_     (1,1)   string  {mustBeMember(fc_, ["r","g","b"])}
                comps_  (1,:)   double  {mustBeNonnegative, mustBeInteger}
                withcm_ (1,1)   logical = true
            end

            this.opts = opts_;
            this.cidx = find(fc_ == this.image_src.MetaData.cOrder);

            % estimate the raw fluorescence(remove background)
            [~, floc] = ismember(comps_, this.mask.MaskComps);
            this.fl(floc, :) = ...
                estimateFluorescence(this.image_src, this.mask.MaskVol, opts_, comps_, fc_);

            if ~opts_.Options.FluorescenceOnly
                % estimate the baseline
                this.bl = estimateBaseline(this.fl, this.opts);

                % calculate the delta F over F
                % omit the neuropil status in 1-photon calcium imaging
                this.dff = (this.fl - this.bl)./this.bl;

                % estimate the noise
                [this.dff, this.noise] = estimateNoise(this.dff, this.opts);
            end

            if withcm_ == true
                this.cmcalc(floc);
            end
        end

        function save(this, comps, folder, fname, withcm)
            % This function saves the data to the given folder
            [~, floc] = ismember(comps, this.mask.MaskComps);
            F = this.Activities(floc, :)';

            % save the neurons response as .xlsx
            T = array2table([this.image_src.Time, F], ...
                "VariableNames", ["Time", this.mask.MaskLabel(floc)]);

            [~, file, ~] = fileparts(fname);
            file_sheet = file + "_activity.xlsx";

            writetable(T, fullfile(folder, file_sheet), "FileType","spreadsheet", ...
                "WriteMode","replacefile");

            if withcm == true
                % save caiman data if possible
                file_cm = fullfile(folder, file+"_cm.mat");
                A = this.A_caiman;
                C = this.C_caiman;
                b = this.b_caiman;
                f = this.f_caiman;
                save(file_cm, "A", "C", "b", "f", '-v7');
            end
        end

        function r = get.Activities(this)
            if this.opts.Options.FluorescenceOnly
                r = this.fl;
            else
                r = this.dff;
            end
        end

        function r = get.Baseline(this)
            r = this.bl;
        end

        function r = get.Noise(this)
            r = this.noise;
        end

        function r = get.CaImAnData(this)
            r = struct("A", this.A_caiman, ...
                       "C", this.C_caiman, ...
                       "b", this.b_caiman, ...
                       "f", this.f_caiman);
        end

        function r = get.Significant(this)
            % use significant level calculate the value
            if ~this.opts.Options.FluorescenceOnly
                switch this.opts.Options.NoiseModel
                    case "normal"
                        switch this.opts.BaselineModel
                            case "MovingQuantile"
                                npts_start = 1;
                            case "MixedExponential"
                                % there may exist bad estimation on the front
                                % of dff because of experiment configurations
                                npts_start = round(size(this.dff, 2)*0.1);
                            otherwise
                        end
                        noi = this.noise(:, npts_start:end);
                        df99 = quantile(this.dff(:, npts_start:end), 0.99, 2);
                        psnr_th = str2double(this.opts.Options.PSNR_Threshold);
                        psnr = 20*log10(abs(df99)./std(noi, 0, 2));
                        r = (psnr > psnr_th);
                    otherwise
                        % no noise model, 100% pass
                        r = true(size(this.fl, 1), 1);
                end
            else
                r = false(size(this.fl, 1), 1);
            end
        end

    end

    methods(Access=private)
        function cmcalc(this, sltidx)
            % This function calculate caiman initial data
            % basic formula: Y = AC + B + E
            % where B can be rewritten as B = bf, where b, f are low-rank
            % matrix when 2p imaging, in this case, for initialization, we
            % treat it as rank one matrix even 1p imaging

            % A comes from mask as initialization
            % flatten the mask
            this.A_caiman = regmask.flatten(this.mask);
            this.A_caiman = this.A_caiman(:, sltidx);
            % renormalization A
            this.A_caiman = this.A_caiman./sum(this.A_caiman, 1);

            % set f as rank one vector with only 1
            this.f_caiman = repmat(1e4, 1, this.image_src.MetaData.frames);
            
            % set b format sparse matrix with rank one
            this.b_caiman = repmat(1e-4, prod(this.mask.size()), 1);
            this.b_caiman = sparse(this.b_caiman);
           
            % calculate the C matrix with CNMF formation
            % use parallel pool for speed up
            src = this.image_src;
            fc = this.cidx;
            Q = (this.A_caiman'*this.A_caiman)\this.A_caiman';
            cbkg = this.opts.Options.Background;    % camera fixed background estimation
            c_caiman = zeros(size(this.A_caiman, 2), this.image_src.MetaData.frames);

            Estimator.auto_parpool("on");

            parfor t = 1:this.image_src.MetaData.frames
                K = double(src.Movie(:,:,fc,:,t) - cbkg);   %#ok<PFBNS>
                K = reshape(K, [], 1);  % flatten 
                c_caiman(:, t) = full(Q) * K;
            end

            Estimator.auto_parpool("off");

            this.C_caiman = c_caiman;
        end
    end

    methods(Static)

        function parobj = auto_parpool(state, nw)
            % This function auto configure the parpool,'thread' has high
            % priority (>=R2022b)
            arguments
                state   (1,1)   string  {mustBeMember(state, ["on","off"])}
                nw      (1,1)   double  {mustBePositive, mustBeInteger} = 8
            end

            switch state
                case "on"
                    % configure if parpool exists
                    parobj = gcp("nocreate");
                    if ~isempty(parobj)
                        if nw == parobj.NumWorkers
                            % keep state
                            return;
                        else
                            delete(parobj);
                        end
                    end

                    if isMATLABReleaseOlderThan("R2022b")
                        parobj = parpool("local", nw);       % older <Processes>
                    else
                        parobj = parpool("Threads", nw);     % regionprops new support
                    end
                case "off"
                    delete(gcp("nocreate"));
                otherwise
            end
        end
    end
end
