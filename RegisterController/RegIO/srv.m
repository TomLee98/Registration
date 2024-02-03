function srv(mov_, v_, mode_, color_)
%SRV This function set reference volumes or raw volumes to mov_
% Input:
%   - mov_: 1-by-1 regmov object
%   - v: m-by-n-by-p(-by-q) volume(s), dimension order must be X,Y,Z(,T)
%   - mode_: 1-by-2 string array, format as [method, frames], where method
%           could be "min","max","mean","median" or "none", etc. Note that
%           "none" with keep the cropped volumes
%   - color_: 1-by-1 string, indicating extract color channel, must be member
%           in ["r", "g", "b"]
%   - standard_: 1-by-1 logical, indicating whether input v follow the
%           standard dimension order: X,Y,Z(,T) or not, true as default
% Output:
%   none
% Behaviour:
%   1. if v_ is a 3D volume, mov_ will be filled with v_ at mode_(2)
%   2. if v_ is 4D volume series, 
%     2.1 if mode_(1) is member in
%     ["min","max","mean","median"], v_ will be processed as v' at first, 
%     then mov_ will be filled with v'. 
%     2.2 if mode_(1) is "none", mov_ will be filled with v_, 
%     in this case, frames of v_ must be equal to mode_(2)

arguments
    mov_        (1,1)       regmov
    v_          (:,:,:,:)
    mode_       (1,2)       string
    color_      (1,1)       string  {mustBeMember(color_, ["r","g","b"])}
end

mf = str2func(mode_(1));
t_range = "[" + string(mode_(2))+ "]";
c_range = string(find(mov_.MetaData.cOrder == color_));
t_loc = (mov_.MetaData.dimOrder=="T");
c_loc = (mov_.MetaData.dimOrder=="C");

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

if ndims(v_) == 3
    % nothing
elseif ndims(v_) == 4
    switch mode_(1)
        case {'min','max'}
            v_ = mf(v_, [], t_loc); %#ok<*NASGU>
        case {'mean','median'}
            v_ = mf(v_, t_loc);
        otherwise
    end
else
    throw(MException("srv:invalidDimensionNumber", ...
        "Unsupported dimensions data."))
end

expr = "mov_.Movie(" + expr + ") = v_;";
eval(expr);
end