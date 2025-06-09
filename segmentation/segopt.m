classdef segopt < handle
    %SEGOPT This class defined as the segmentor options type

    properties(Access = public)
        %% common
        method      (1,1)   string  {mustBeMember(method, ["morphology", "cellpose"])} = "morphology"
        modal       (1,1)   string  {mustBeMember(modal, ["train", "predict"])} = "predict"
        bkg         (1,1)   double  {mustBeNonnegative}                 = 100
        gamma       (1,1)   double  {mustBeInRange(gamma, 0, 2)}        = 1.0
        bdbox       (:,6)   double  {mustBeNonnegative, mustBeInteger}  = []

        %% morphology
        r_u         (1,1)   double  {mustBePositive}                    = 2.2
        r_s         (1,1)   double  {mustBePositive}                    = 0.6
        n_max       (:,1)   double  {mustBePositive, mustBeInteger}     = 100
        fg_th       (1,1)   double  {mustBeInRange(fg_th, 0, 255)}      = 127
        auto_fg     (1,1)   logical                                     = true
        sph_th      (1,1)   double  {mustBeInRange(sph_th, 0, 1)}       = 0.5
        rm_outlier  (1,1)   string  {mustBeMember(rm_outlier, ["left","right","both"])} = "both"
        ws_old      (1,1)   logical                                     = false

        %% machine learning (cellpose)
        r_mean      (1,1)   double  {mustBePositive}                    = 1.8
        r_min       (1,1)   double  {mustBePositive}                    = 1.4
        model       (1,1)   string  {mustBeMember(model,["cyto","CP","nuclei","KC"])} = "nuclei"
        stage       (1,1)   string  {mustBeMember(stage, ["L1", "L2", "L3"])} = "L1"
        cell_th     (1,1)   double  {mustBeInRange(cell_th, -6, 6)}     = -3
        enhance_ml  (1,1)   logical                                     = false
        batch_size  (1,1)   double  {mustBeInRange(batch_size, 1, 8),  mustBeInteger} = 4
        learn_rate  (1,1)   double  {mustBeInRange(learn_rate, 0, 1)}   = 0.01
        weight_decay(1,1)   double  {mustBeInRange(weight_decay, 0, 1)} = 1E-5
        epochs_max  (1,1)   double  {mustBePositive, mustBeInteger}     = 500
        acceleration(1,1)   string  {mustBeMember(acceleration, ["openvino", "none"])} = "openvino"
    end
    
    methods
        function this = segopt()
            %SEGOPT A constructor
        end

    end
end

