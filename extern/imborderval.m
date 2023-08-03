function v = imborderval(A, varargin)
%GENFILLVAL This function generate border value for 2D or 3D image
% input:
%   - A: 2-D or 3D image
%   - varargin:
%       - (optional) r: the border size vector, if A is 2D image, r could 
%                       be vector with 4 elements, [left, bottom, right, top],
%                       if A is 3D image, r is vector with 6 elements, 
%                       [left, bottom, right, top, front, back], which also
%                       could be scalar, this will apply on each direction,
%                       2 as default
%       - (optional) method: the border calculate method at each direction, 
%                            which could be "min", "max", "mean", "median",
%                            "other", "mean" as default
%       - quantile: the Name-Value pair for quantile calculation, which
%                   need a quantile number p
%       - selfdef: the Name-Value pair for self defination funcition
% output:
%   - v: the border vector which has the same size with r

VALID_METHOD = ["min","max","mean","median","other"];

p = inputParser();
valid_A = @(x)validateattributes(x,{'numeric'},{'real'});
valid_r = @(x)validateattributes(x,{'numeric'},{'vector'});
valid_method = @(x)ismember(x, VALID_METHOD);
valid_quantile = @(x) validateattributes(x,'numeric',{'scalar','>=',0,'<=',1});
valid_selfdef = @(x) validateattributes(x,{'function_handle'},{'scalar'});

default_r = 2;
default_method = "mean";
default_quantile = 0;
default_selfdef = @function_handle.empty;

addRequired(p,"A",valid_A);
addOptional(p,"r",default_r,valid_r);
addOptional(p,"method",default_method,valid_method);
addParameter(p,"quantile",default_quantile,valid_quantile);
addParameter(p,"selfdef",default_selfdef,valid_selfdef);
parse(p, A, varargin{:});

size_A = size(A);

if ismatrix(p.Results.A) || isrgb(p.Results.A)
    % 2D problem
    assert(numel(p.Results.r)==1 || numel(p.Results.r)==4, ...
        "imborderval:borderCountError","Image border number is invalid.");
    if numel(p.Results.r) == 1
        r = repmat(p.Results.r, 1, 4);
    else
        r = p.Results.r;
    end

    if r(1)+r(3)+1 >= size_A(2) || r(5)+r(6)+1 >= size_A(1)
        error("imborderval:borderOverflow","Image border is too large")
    end
elseif ndims(p.Results.A) == 3
    % 3D problem
    assert(numel(p.Results.r)==1 || numel(p.Results.r)==6, ...
        "imborderval:borderCountError","Image border number is invalid.");
    if numel(p.Results.r) == 1
        r = repmat(p.Results.r, 1, 6);
    else
        r = p.Results.r;
    end

    if r(1)+r(3)+1 >= size_A(2) || r(5)+r(6)+1 >= size_A(1) ||...
            r(2)+r(4)+1 >= size_A(3)
        error("imborderval:borderOverflow","Image border is too large")
    end
end

% calculate the border value
v = zeros(size(r));

% cut the border

if numel(r) == 4
    border_left = A(:, 1:r(1));
    border_bottom = A(end-r(2)+1:end, r(1)+1:end-r(3));
    border_right = A(:, end-r(3)+1:end);
    border_top = A(1:r(4), r(1)+1:end-r(3));
    border_front = [];
    border_back = [];
elseif numel(r) == 6
    border_left = A(:, 1:r(1), :);
    border_bottom = A(:, r(1)+1:end-r(1), 1:r(2));
    border_right = A(:, end-r(3)+1:end, :);
    border_top = A(:, r(1)+1:end-r(1), end-r(4)+1:end);
    border_front = A(end-r(5)+1:end, r(1)+1:end-r(3), r(2)+1:end-r(4));
    border_back = A(1:r(6), r(1)+1:end-r(3), r(2)+1:end-r(4));
end

switch p.Results.method
    case "other"
        if p.Results.quantile ~= 0
            v(1) = quantile(border_left, p.Results.quantile, "all");
            v(2) = quantile(border_bottom, p.Results.quantile, "all");
            v(3) = quantile(border_right, p.Results.quantile, "all");
            v(4) = quantile(border_top, p.Results.quantile, "all");
            v(5) = quantile(border_front, p.Results.quantile, "all");
            v(6) = quantile(border_back, p.Results.quantile, "all");
            return;
        else
            if p.Results.selfdef ~= @function_handle.empty
                f = p.Results.selfdef;
            end
        end
    otherwise
        f = str2func(p.Results.method);
end

% calculate the border value
v(1) = f(border_left,"all");
v(2) = f(border_bottom,"all");
v(3) = f(border_right,"all");
v(4) = f(border_top,"all");
v(5) = f(border_front,"all");
v(6) = f(border_back,"all");
end


function r = isrgb(A)
if isnumeric(A)
    if ndims(A) == 3 && size(A,3) == 3
        switch class(A)
            case 'double'
                if all(A > 0 & A < 1)
                    r = true;
                else
                    r = false;
                end
            otherwise
                r = true;
        end
    else
        r = false;
    end
else
    r = false;
end
end