classdef Estimator < handle
    %SIGNALEXTRACTOR This class define a signal extractor object

    properties(Constant, Hidden)

    end

    properties(GetAccess=public, Dependent)
        Activities      % get,  variable

    end
    
    properties(Access=private, Hidden)
        activities      % n-by-t double array, row for components, colume for time step
        mask            % m-by-n-by-p nonnegtive integer label array
        lbl_mask        % 1-by-s positive integer, labels in label array
    end

    properties(SetAccess=immutable, Hidden)
        image_src      % 1-by-1 regmov object, the calcium image data 
    end
    
    methods
        function this = Estimator(src_, mask_)
            arguments
                src_     (1,1)   regmov
                mask_    (:,:,:) 
            end

            this.image_src = src_;
            this.mask = mask_;
            this.lbl_mask = unique(this.mask);
            this.lbl_mask(this.lbl_mask==0) = [];
            ncomp = numel(this.lbl_mask);

            this.activities = nan(ncomp, this.image_src.MetaData.frames);
        end
        
        function fit(this, opts_, ch_, comps_)
            % This function fit the signal model
            arguments
                this
                opts_   (1,1)   sigopt
                ch_     (1,1)   string  {mustBeMember(ch_, ["r","g","b"])}
                comps_  (1,:)   double  {mustBeNonnegative, mustBeInteger}
            end

            fn = this.image_src.MetaData.frames;
            c = (ch_ == this.image_src.MetaData.cOrder);
            acts = nan(numel(comps_), fn);
            src_img = this.image_src;
            mask_ = this.mask;

            % use parallel for speed up
            delete(gcp("nocreate"));
            parpool("Threads");

            %TODO: use util estimators

            parfor n = 1:numel(comps_)
                mk = (mask_==comps_(n));

                for t = 1:fn
                    tmpvol = src_img.Movie(:,:,c,:,t); %#ok<PFBNS>
                    acts(n, t) = mean(tmpvol(mk));
                end

            end

            [~, locs] = ismember(comps_, this.lbl_mask);
            this.activities(locs, :) = acts;

            % close parpool
            delete(gcp("nocreate"));
        end

        function r = get.Activities(this)
            r = this.activities;
        end

    end
end

