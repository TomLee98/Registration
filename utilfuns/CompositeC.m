function C = CompositeC(R, G, B, cdef)
%COMPOSITEC_NEW composites the color channel and generates the RGB
% representation
% Input:
%   - R: 2-D/3-D uint8 matrix for channel R
%   - G: 2-D/3-D uint8 matrix for channel G
%   - B: 2-D/3-D uint8 matrix for channel B
%   - cdef: 3-by-3 color transformation array, 0 ~ 1
% Output:
%   - C: 3-D/4-D uint8 true color matrix

% Version 1.0.0
% Copyright (c) 2022-2025, Weihan Li

arguments
    R       (:,:,:) single   {isnumeric(R)}
    G       (:,:,:) single   {isnumeric(G)}
    B       (:,:,:) single   {isnumeric(B)}
    cdef    (3,3)   double  {mustBeInRange(cdef, 0, 1)}
end

% validate image size
if isempty(R) && isempty(G) && isempty(B)
    C = uint8([]);
    return;
end
if (~isempty(R)&&~isempty(G)&&(any(size(R,1:3)~=size(G,1:3)))) ...
        || (~isempty(G)&&~isempty(B)&&any(size(G,1:3)~=size(B,1:3)))
    throw(MException("Composite:sizeNotMatch", ...
        "Input channels size are not compatible."));
end

% compatible modification
vsize = max(max(size(R, 1:3), size(G, 1:3)), size(B, 1:3));
vblack = zeros(vsize, "single");
if isempty(R), R = vblack; end
if isempty(G), G = vblack; end
if isempty(B), B = vblack; end

% transform 
Rp = cdef(1,1)*R + cdef(1,2)*G + cdef(1,3)*B;
Gp = cdef(2,1)*R + cdef(2,2)*G + cdef(2,3)*B;
Bp = cdef(3,1)*R + cdef(3,2)*G + cdef(3,3)*B;

% explicit casting to uint8
C = cast(cat(ndims(vblack)+1, Rp, Gp, Bp), "uint8");
end

