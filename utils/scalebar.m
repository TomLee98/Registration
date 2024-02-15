classdef scalebar < handle
    %SCALEBAR This class define a scale bar object, which could push in
    % uifigure containter

    properties(Constant, Hidden)
        DRAW_SUCCESS = 0
        UPDATE_SUCCESS = 0
        DRAW_FAILED = -1
    end
    
    properties(Access=private, Hidden)
        handle_x
        handle_y
        txtobj_x
        txtobj_y
        scale
        sbar
        font
        ismoving
        islocked
        sf
        imwidth
        imheight
    end

    properties(Access=public, Dependent)
        Scale
        Bar
        Font
        ScaleFactor
        Locked
        Moving
    end

    properties(SetAccess=immutable, GetAccess=public)
        Parent
    end
    
    methods
        function this = scalebar(parent, w, h)
            arguments
                parent  (1,1)   matlab.ui.control.UIAxes
                w       (1,1)   double  {mustBePositive, mustBeInteger} = 1024
                h       (1,1)   double  {mustBePositive, mustBeInteger} = 1024
            end
            %SCALEBAR A Constructor
            this.scale = struct("xr",[], ...
                                "yr",[], ...
                                "unit","um");
            this.sbar = struct("xlen",10, ...
                               "ylen",10, ...
                               "ishor",false, ...
                               "isver",false, ...
                               "linw",1, ...
                               "color","w", ... 
                               "mk","none", ...
                               "sx",40, ...
                               "sy",40,...
                               "mk_fill",true, ...
                               "mk_size",3, ...
                               "txt_space",8);
            this.font = struct("FontName",'Arial', ...
                               "FontWeight",'bold', ...
                               "FontAngle",'normal', ...
                               "FontUnits",'points', ...
                               "FontSize",10);

            this.ismoving = false;
            this.islocked = false;
            this.sf = 1;
            this.imwidth = w;
            this.imheight = h;

            this.Parent = parent;
        end
        
        function r = get.Scale(this)
            r = this.scale;
        end

        function set.Scale(this, r_)
            arguments
                this
                r_  (1,1)   struct  {mustBeScale}
            end

            this.scale = r_;
        end

        function r = get.Bar(this)
            r = this.sbar;
        end

        function set.Bar(this, r_)
            arguments
                this
                r_  (1,1)   struct  {mustBeBar}
            end

            this.sbar = r_;
        end

        function r = get.Font(this)
            r  = this.font;
        end

        function set.Font(this, r_)
            arguments
                this
                r_      (1,1)   struct  {mustBeFont}
            end

            this.font = r_;
        end

        function r = get.ScaleFactor(this)
            r = this.sf;
        end

        function set.ScaleFactor(this, r_)
            arguments
                this
                r_  (1,1)   double  {mustBePositive}
            end

            this.sf = r_;
        end

        function r = get.Locked(this)
            r = this.islocked;
        end

        function set.Locked(this, r_)
            arguments
                this
                r_  (1,1)   logical
            end

            this.islocked = r_;
        end

        function r = get.Moving(this)
            r = this.ismoving;
        end

        function set.Moving(this, r_)
            arguments
                this
                r_  (1,1)   logical
            end

            this.ismoving = r_;
        end
    end

    methods(Access=public)
        function status = draw(this)
            if ((this.sbar.xlen/this.scale.xr >= this.imwidth - 20) && this.sbar.ishor)...
                    || ((this.sbar.ylen/this.scale.yr >= this.imheight - 20) && this.sbar.isver)
                status = this.DRAW_FAILED;
                return;
            end

            % update the selected bar
            if this.sbar.mk_fill == true
                fc = this.sbar.color;
            else
                fc = "none";
            end

            if this.sbar.ishor == true
                switch this.sbar.mk
                    case "flat"
                        ept = "|";
                    otherwise
                        ept = this.sbar.mk;
                end

                % generate the horizontal scale bar
                this.handle_x = line(this.Parent, ...
                    [this.sbar.sx, this.sbar.sx+this.sbar.xlen/this.scale.xr], ...
                    [this.sbar.sy, this.sbar.sy], ...
                    "Color", this.sbar.color,"LineWidth",this.sbar.linw,"Marker",ept,...
                    "MarkerSize",this.sbar.mk_size,"MarkerFaceColor",fc);
                % binding text on scale bar
                this.txtobj_x = text(this.Parent,...
                    this.sbar.sx+this.sbar.xlen/this.scale.xr/2,...
                    this.sbar.sy+this.sf*this.sbar.txt_space,...
                    string(this.sbar.xlen)+" "+this.scale.unit,'Color',this.sbar.color,'FontSize',this.font.FontSize,...
                    'FontWeight',this.font.FontWeight,'FontName',this.font.FontName,...
                    'FontAngle',this.font.FontAngle,'FontUnits',this.font.FontUnits,...
                    'HorizontalAlignment','center');
                % add self define callback
                this.handle_x.ButtonDownFcn = ...
                    @(src, evt)this.ButtonDownFunc(src, evt);

                status = this.DRAW_SUCCESS;
            end

            if this.sbar.isver == true
                switch this.sbar.mk
                    case "flat"
                        ept = "_";
                    otherwise
                        ept = this.sbar.mk;
                end
                % generate the vertical scale bar
                this.handle_y = line(this.Parent, ...
                    [this.sbar.sx, this.sbar.sx]+this.sbar.xlen/this.scale.xr, ...
                    [this.sbar.sy, this.sbar.sy-this.sbar.ylen/this.scale.yr], ...
                    "Color", this.sbar.color,"LineWidth",this.sbar.linw,"Marker",ept,...
                    "MarkerSize",this.sbar.mk_size,"MarkerFaceColor",fc);
                % binding text on scale bar
                this.txtobj_y = text(this.Parent,...
                    this.sbar.sx+this.sbar.xlen/this.scale.xr+this.sf*this.sbar.txt_space,...
                    this.sbar.sy-this.sbar.ylen/this.scale.yr/2,...
                    string(this.sbar.ylen)+" "+this.scale.unit,'Color',this.sbar.color,'FontSize',this.font.FontSize,...
                    'FontWeight',this.font.FontWeight,'FontName',this.font.FontName,...
                    'FontAngle',this.font.FontAngle,'FontUnits',this.font.FontUnits,...
                    'Rotation',90,'HorizontalAlignment','center');
                % add self define callback
                this.handle_y.ButtonDownFcn = ...
                    @(src, evt)this.ButtonDownFunc(src, evt);
            end
        end

        function update(this, metadata)
            arguments
                this
                metadata    (1,12)  table
            end

            this.imwidth = metadata.width;
            this.imheight = metadata.height;
            this.scale.xr = metadata.xRes;
            this.scale.yr = metadata.yRes;
            this.sbar.xlen = min(this.sbar.xlen, 0.5*metadata.width*metadata.xRes);
            this.sbar.ylen = min(this.sbar.ylen, 0.5*metadata.height*metadata.yRes);
            this.sbar.sy = metadata.height - this.sbar.sx;
            this.sf = min(get(0).ScreenSize(3:4)./[metadata.width, metadata.height]);
        end

        function move(this, mpcur)
            arguments
                this
                mpcur   (1,2)   double
            end

            bar_len_x = this.sbar.xlen/this.scale.xr;
            bar_len_y = this.sbar.ylen/this.scale.yr;
            if ~isempty(this.handle_x) && isvalid(this.handle_x)
                this.handle_x.XData = [mpcur(1)-bar_len_x/2, ...
                                       mpcur(1)+bar_len_x/2];
                this.handle_x.YData = [mpcur(2), mpcur(2)];
                this.txtobj_x.Position = [mpcur(1), ...
                    mpcur(2)+this.sf*this.sbar.txt_space, 0];
            end
            if ~isempty(this.handle_y) && isvalid(this.handle_y)
                this.handle_y.XData = [mpcur(1)+bar_len_x/2, ...
                                       mpcur(1)+bar_len_x/2];
                this.handle_y.YData = [mpcur(2), mpcur(2)-bar_len_y];
                this.txtobj_y.Position = [mpcur(1)+bar_len_x/2 ...
                    + this.sf*this.sbar.txt_space, mpcur(2)-bar_len_y/2, 0];
            end
        end

        function delete(this)
            rmobjs(this);

            % ~
        end

        function clear(this)
            % remove objects on screen
            rmobjs(this);
        end
    end

    methods(Access=private, Hidden)
        function ButtonDownFunc(this, src, event)
            src; %#ok<VUNUS>
            % This function is button down callback
            if this.islocked == true
                return;
            end
            switch event.Button
                case 1
                    % left button for motion tracking
                    % switch the moving status
                    this.ismoving = ~this.ismoving;
                    % lazy update sx, sy
                    if this.ismoving == false
                        mp = this.Parent.CurrentPoint(1,1:2);
                        this.sbar.sx = mp(1);
                        this.sbar.sy = mp(2);
                    end
                otherwise
            end
        end

        function rmobjs(this)
            % remove the exist scale bar
            if ~isempty(this.handle_x) && isvalid(this.handle_x)
                delete(this.handle_x);
                this.handle_x = [];
            end
            if ~isempty(this.handle_y) && isvalid(this.handle_y)
                delete(this.handle_y);
                this.handle_y = [];
            end
            if ~isempty(this.txtobj_x) && isvalid(this.txtobj_x)
                delete(this.txtobj_x);
                this.txtobj_x = [];
            end
            if ~isempty(this.txtobj_y) && isvalid(this.txtobj_y)
                delete(this.txtobj_y);
                this.txtobj_y = [];
            end
        end
    end
end

function mustBeScale(A)
fields = fieldnames(A);
if ~isempty(setxor(fields, ["xr","yr","unit"]))
    throw(MException("mustBeScale:invalidArguments", ...
        "Invalid scale structure."));
end

validateattributes(A.xr, "double", {'positive','scalar'});
validateattributes(A.yr, "double", {'positive','scalar'});
validateattributes(A.unit, "string", {'scalar'});
end

function mustBeBar(A)
fields = fieldnames(A);
if ~isempty(setxor(fields, ["xlen","ylen","ishor","isver","linw","color", ...
        "mk", "sx", "sy", "mk_fill", "mk_size", "txt_space"]))
    throw(MException("mustBeBar:invalidArguments", ...
        "Invalid bar structure."));
end

% omit value validation
end

function mustBeFont(A)
fields = fieldnames(A);

if ~isempty(setxor(fields, ["FontName","FontWeight","FontAngle", ...
        "FontUnits","FontSize"]))
    throw(MException("mustBeFont:invalidArguments", ...
        "Invalid font structure."));
end

end