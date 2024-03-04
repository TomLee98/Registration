classdef regtmpl
    %REGTMPL This class is registration template defination, which supports

    properties(SetAccess=immutable, GetAccess=private, Hidden)
        regdata
        fixvdef
        auto_update
    end

    properties(Access=private, Hidden)
        refvol
    end

    properties(GetAccess=public, Dependent)
        RefVol
    end
    
    methods
        function this = regtmpl(regmov_, fixvdef_, auto_update_)
            %REGTMPL A Constructor
            arguments
                regmov_     (1,1) regmov {mustBeTemplatable}
                fixvdef_    (1,1) struct {mustBeFixedVolDefination} % with 'Global','Local','Channel'
                auto_update_(1,1) logical = false
            end
            
            this.regdata = regmov_;     % shallow copy handle
            this.fixvdef = fixvdef_;
            this.auto_update = auto_update_;    % perfermance option

            this = update_ref(this);
        end
        
        function r = get.RefVol(this)
            % case off for fast reading without update
            if this.auto_update
                % update the reference volume
                this = update_ref(this);
            end

            r = this.refvol;
        end
    end

    methods(Access=private, Hidden)
        function this = update_ref(this)
            % Note that: dimension order standarded as XYZ
            this.refvol = grv(this.regdata, ...
                              this.fixvdef.Sampling, ...
                              this.fixvdef.Channel);
        end
    end
end


% ====================== local utility function ===================
function mustBeFixedVolDefination(A)
VALID_METHOD = ["min","max","mean","median"];
fields = fieldnames(A);
if ~isempty(setxor(fields, ["Sampling", "Channel"]))
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Template defination must be with 'Sampling' and 'Channel' field."));
end
if (numel(A.Sampling) ~= 2) || ~isstring(A.Sampling)
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "'Sampling' must be with two string elements: <method>, <frames>."));
end
if ~all(ismember(A.Channel, ["r","g","b"]))
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Channel can only be member in [""r"",""g"",""b""]"));
end
if ~ismember(lower(A.Sampling(1)), VALID_METHOD)
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Only 'min','max','mean','median' are valid template generation methods."));
end
if ~isvector(str2num(A.Sampling(2))) %#ok<ST2NM>
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Frames defination is non-convertible."));
end
end

function mustBeTemplatable(A)
if ~ismember("T", A.MetaData.dimOrder) ...
        || A.MetaData.frames == 0
    throwAsCaller(MException("regmpl:mustBeTemplatable:" + ...
        "invalidRegistrationMovie", ...
        "Can not generate template because time dimension lost."));
end
end