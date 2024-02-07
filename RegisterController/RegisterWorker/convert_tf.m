function B = convert_tf(A, varargin)
%CONVERT_TF This function convert geometric transformation to satisfy the
%imregtform initialization requirement
% Input:
%   - A: a geometric transformation, must be affine3d or affinetform3d
%        object, which can be translated to rigid or translated
%   - type: 1-by-1 string, could be "translation", "rigid" or "affine"
% Output:
%   - B: the translated geometric transformation object, could be "affine3d",
%       "rigid3d", or "affinetform3d", "rigidtform3d" or "transltform3d",
%       which depends on your input formation

p = inputParser();
addRequired(p, "A", @(x)validateattributes(x, {'affine3d', 'affinetform3d'}, "scalar"));
addOptional(p, "type_", "affine", @(x)validateattributes(x, {'string', 'char'}, "scalartext"));
p.parse(A, varargin{:});

A = p.Results.A;
type_ = p.Results.type_;
if ~ismember(type_, ["affine","rigid", "translation"])
    throw(MException("convert_tf:invalidConvertType", ...
        "Undefined transformation geometric type."));
end

switch type_
    case "translation"
        if isa(A, "affine3d")
            if isTranslation(A)
                B = A;
            else
                % omit the scale change
                A.T(1:A.Dimensionality, 1:A.Dimensionality) ...
                         = eye(A.Dimensionality);
                B = A;
            end
        else
            A.T(1:A.Dimensionality, 1:A.Dimensionality) ...
                = eye(A.Dimensionality);
            B = A;
        end
    case "rigid"
        if isa(A, "affine3d")
            
        else

        end
    case "affine"
        B = A;
    otherwise
end

end

