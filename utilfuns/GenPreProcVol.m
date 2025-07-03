function v = GenPreProcVol(vol, opt_reg)
%GENPREPROCVOL This function generate pre-processed volume by given raw
% 3-D array and registration option
arguments
    vol         (:,:,:)     uint16
    opt_reg     (1,1)       regopt
end

switch opt_reg.Algorithm
    case "OCREG"
        switch opt_reg.Mode
            case "global"
                ds = opt_reg.Options.DS;
                if ds == "auto", regds = inf;
                else,            regds = 1/str2double(ds.extractBefore("X"));
                end
                [~, vol] = Resample(vol, regds);
                v = preproc_oc(vol, ...
                    opt_reg.Options.MedianFilter, ...
                    opt_reg.Options.GaussianFilter);
                
                % v = isoutlier(single(v), opt_reg.Options.CoarseArgs.Filter);
                % v = bwareaopen(v, opt_reg.Options.CoarseArgs.VT, 26);
                % v = imerode(v, strel("sphere", opt_reg.Options.CoarseArgs.Radius));
            otherwise
                v = vol;
        end
    case "TCREG"
        switch opt_reg.Mode
            case "global"
                ds = opt_reg.Options.DS;
                if ds == "auto", regds = inf;
                else,            regds = 1/str2double(ds.extractBefore("X"));
                end
                [~, vol] = Resample(vol, regds);
                v = preproc_tc(vol, ...
                    opt_reg.Options.DilateFilter, ...
                    opt_reg.Options.MedianFilter, ...
                    opt_reg.Options.GaussianFilter, ...
                    opt_reg.Options.Gamma);
            otherwise
                v = vol;
        end
    case "LTREG"
        ds = opt_reg.Options.DS;
        if ds == "auto", regds = inf;
        else,            regds = 1/str2double(ds.extractBefore("X"));
        end
        [~, vol] = Resample(vol, regds);
        v = preproc_tc(vol, opt_reg.Options.DilateFilter);
    otherwise
        v = vol;
end

end

