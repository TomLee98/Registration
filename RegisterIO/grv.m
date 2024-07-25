function v = grv(mov_, mode_, color_, standard_)
% GRV   This function generate reference volume or extract some frames
%without processing at color_ channel
% Input:
%   - mov_: 1-by-1 regmov object
%   - mode_: 1-by-2 string array, format as [method, frames], where method
%           could be "min","max","mean","median" or "none", etc. Note that
%           "none" with keep the cropped volumes
%   - color_: 1-by-1 string, indicating extract color channel, must be member
%           in ["r", "g", "b"]
%   - standard_: 1-by-1 logical, indicating whether output v follow the
%           standard dimension order: X,Y,Z(,T) or not, true as default
% Output:
%   - v: m-by-n-by-p(-by-t) volume(s)

arguments
    mov_        (1,1)   regmov
    mode_       (1,2)   string
    color_      (1,1)   string  {mustBeMember(color_, ["r","g","b"])}
    standard_   (1,1)   logical = true
end

mf = str2func(mode_(1));
t_range = "[" + string(mode_(2))+ "]";
c_range = string(find(mov_.MetaData.cOrder == color_));
t_loc = find(mov_.MetaData.dimOrder=="T");
c_loc = find(mov_.MetaData.dimOrder=="C");

if isempty(c_range) || isempty(t_loc) || isempty(c_loc)
    throw(MException("grv:invalidMovie", ...
        "Bad calling: invalid movie dimension."));
end

mov_ndim = numel(mov_.MetaData.dimOrder);
expr = "";
for dp = 1:mov_ndim
    if dp == c_loc
        expr = expr + c_range;
    elseif dp == t_loc
        expr = expr + t_range;
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
        v = D;
end
v = squeeze(v);

if ismatrix(v) || ndims(v)> 4
    throw(MException("regtmpl:grv:innerError", ...
        "Invalid reference volume generated."));
end

if standard_ == true
    % "C" is selected by color_
    [~, p] = ismember(["X","Y","Z","T"], mov_.MetaData.dimOrder);
    [~, vp] = sort(p, "ascend");
    v = permute(v, vp);
end

% change the value type to uint16
v = uint16(v);
end

