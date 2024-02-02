classdef regtmpl
    %REGTMPL This class is registration template defination, which supports
    
    properties(SetAccess=immutable, GetAccess=private, Hidden)
        regdata
        fixvdef
    end

    properties(Access=private, Hidden)
        refvol
    end

    properties(GetAccess=?RegisterWorker, Dependent)
        RefVol
    end
    
    methods
        function this = regtmpl(regmov_, fixvdef_)
            %REGTMPL A Constructor
            arguments
                regmov_     (1,1) regmov {mustBeTemplatable}
                fixvdef_    (1,1) struct {mustBeFixedVolDefination} % with 'Global','Local','Channel'
            end
            
            this.regdata = regmov_;     % shallow copy handle
            this.fixvdef = fixvdef_;

            this = gen_refvol(this);
        end
        
        function r = get.RefVol(this)
            r = this.refvol;
        end
    end

    methods(Access=private, Hidden)
        function this = gen_refvol(this)
            this.refvol = struct("Global", [], ...
                                 "Local", []);

            this.refvol.Global = grv(this.regdata, ...
                                     this.fixvdef.Global, ...
                                     this.fixvdef.Channel);
            this.refvol.Local = grv(this.regdata, ...
                                    this.fixvdef.Local, ...
                                    this.fixvdef.Channel);

            function v = grv(mov_, mode_, color_)
                mf = str2func(mode_(1));
                t_range = "[" + string(mode_(2))+ "]";
                c_range = string(find(mov_.MetaData.cOrder == color_));
                t_loc = (mov_.MetaData.dimOrder=="T");
                c_loc = (mov_.MetaData.dimOrder=="C");

                mov_ndim = numel(mov_.MetaData.dimOrder);
                expr = "";
                for dp = 1:mov_ndim
                    if dp == c_loc
                        expr = expr + t_range;
                    elseif dp == t_loc
                        expr = expr + c_range;
                    else
                        expr = expr + ":";
                    end
                    if dp ~= mov_ndim, expr = expr + ","; end
                end
                expr = "mov_.Movie(" + expr + ");";
                D = eval(expr);     % create croped temporary data: D
                switch mode_(1)
                    case {'min','max'}
                        v = mf(D, [], t_loc);
                    case {'mean','median'}
                        v = mf(D, t_loc);
                    otherwise
                end
                v = unique(v);

                if ndims(v) ~= 3
                    throw(MException("regtmpl:grv:innerError", ...
                        "Invalid reference volume generated."));
                end
            end
        end
    end
end


% ====================== local utility function ===================
function mustBeFixedVolDefination(A)
VALID_METHOD = ["min","max","mean","median"];
fields = fieldnames(A);
if ~isempty(setxor(fields, ["Global", "Local", "Channel"]))
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Template defination must be with 'Global', 'Local' and 'Channel' field."));
end
if (numel(A.Global) ~= 2) || (numel(A.Local) ~= 2) ...
        || ~isstring(A.Global) || ~isstring(A.Local)
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "'Global' or 'Local' must be with two string elements: <method>, <frames>."));
end
if ~all(ismember(A.Channel, ["r","g","b"]))
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Channel can only be member in [""r"",""g"",""b""]"));
end
if ~ismember(lower(A.Global(1)), VALID_METHOD) ...
        || ~ismember(lower(A.Local(1)), VALID_METHOD)
    throwAsCaller(MException("regmpl:mustBeFixedVolDefination:" + ...
        "invalidTemplateDefination", ...
        "Only 'min','max','mean','median' are valid template generation methods."));
end
if ~isvector(str2num(A.Global(2))) || ~isvector(str2num(A.Local(2))) %#ok<ST2NM>
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