function C = TransformBCG(A, param)
% this function transform image A (m*n matrix) intensity by
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
    param (1,1) struct = struct("b",0,"c",1,"ga",1);
end

% rescale the pixel intensity
A = rescale(A,0,intmax("uint16"));

% linear model
C = imlincomb(param.c, uint16(single(A).^param.ga), param.b);
C = uint8(rescale(C, 0, intmax("uint8")));
end

