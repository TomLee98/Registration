function C = TransformBCG(A, param)
% this function transform image A (m*n*p matrix) intensity by
% using param struct
% input:
%   - A: the image need to be transformed,
%        double/single/uint8/uint16
%   - param: the tramsform parameters, struct with {b,c,ga},
%            image intensity model: y = c*x^ga+b
% output:
%   - C: the transformed image, uint8 for fast viewing

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A (:,:,:) {mustBeNumeric}   % RGB
    param (1,1) struct = struct("r",[0,65535],"b",0,"c",1,"ga",1);
end

% rectangle constrain
A = max(min(A, param.r(2)), param.r(1));

% rescale the pixel intensity
A = matlab.internal.math.rescaleKernel(A, 0, 65535, param.r(1), param.r(2));

% linear model
C = imlincomb(param.c, uint16(single(A).^param.ga), param.b);
iminc = double(min(C,[],"all")); 
imaxc = double(max(C, [], "all"));
C = cast(matlab.internal.math.rescaleKernel(C, 0, 255, iminc, imaxc), "uint8"); 
end

