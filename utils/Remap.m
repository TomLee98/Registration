function J = Remap(I, type, LUT)
% This function remap the image I intensity to new type by defined LUT and
% it will fully map the target gray scale
arguments
    I {mustBeNonnegative};
    type (1,1) string {mustBeMember(type, ["uint8","uint16","uint32","single","double"])} = "uint8";
    LUT (1,2) double {mustBeNonnegative} = [0, inf];
end

if all(LUT == [0, inf])
    % auto scale the data by min and max
    LUT = [min(I,[],"all"), max(I,[],"all")];
else
    if LUT(1) < min(I,[],"all") || LUT(2) > max(I,[],"all")
        warning("Auto shrink the LUT to fit the gray scale.");
        LUT(1) = max(LUT(1), min(I,[],"all"));
        LUT(2) = min(LUT(2), max(I,[],"all"));
    end
end

validateattributes(LUT,{'numeric'},{'increasing', 'row'});

% remap the intensity
J = double(I - LUT(1))/double(LUT(2)- LUT(1));

if ismember(type, ["uint8","uint16","uint32"])
    J = cast(J*double(intmax(type)), type);
elseif type == "single"
    J = single(J);
end

end

