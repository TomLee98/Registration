function C = TransformBCG(A, param)
% this function transform image A (m*n matrix) intensity by
% using param struct
% input:
%   - A: the image need to be transformed,
%        double/single/uint8/uint16
%   - param: the tramsform parameters, struct with {b,c,ga},
%            image intensity model: y = c*x^ga+b
% output:
%   - C: the transformed image

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A (:,:,:) {mustBeNumeric}   % RGB
    param (1,1) struct = struct("b",0,"c",1,"ga",1);
end

% rescale the pixel intensity
if ismember(class(A),["uint8","uint16","uint32","uint64"])
    A = rescale(A,0,intmax(class(A)));
else
    A = rescale(uint16(A),0,65535); % scale as uint16
end

% linear model
C = imlincomb(param.c, uint16(single(A).^param.ga), param.b);
end

