function A = MaskMovie(A, mask, inverse)
%MASKMOVIE This function masked A by using mask
% Input:
%   - A: the movie need to be masked, 3D/4D/5D array, where 3D array with
%        dimension order(x,y,z), 4D array with (x,y,z,t) and 5D array with
%        (x,y,c,z,t).
%   - mask: the volume mask, must be 3D logical array.
%   - inverse: the flag for mask inverse. By default, true means the pixel
%              is selected and false means hide. If inverse is true, the
%              the mask will be not(mask).
% Output:
%   - A: the masked movie

% Version 1.0.0
% Copyright (c) 2022-2024, Weihan Li

arguments
    A
    mask (:,:,:) logical
    inverse (1,1) logical = false
end

if isempty(A), return; end

D = ndims(A);
assert(D>=3 && D<=5, "MaskMovie:invalidMovieFormat","Invalid dimension of movie.");
if D == 3
    % A dimension order as (x,y,z)
    vol_size_A = size(A);
elseif D == 4
    % A dimension order as (x,y,z,t)
    vol_size_A = size(A,[1,2,3]);
elseif D == 5
    % A dimension order as (x,y,c,z,t)
    vol_size_A = size(A,[1,2,4]);
end
size_mask = size(mask);
assert(all(size_mask==vol_size_A),"MaskMovie:sizeNotMatch",...
    "Mask is not match.");

if inverse, mask = ~ mask; end %#ok<NASGU> 

if class(A) ~= "mpimg"
    mask = eval([class(A),'(mask)']);
else
    mask = eval([A.DataType, '(mask)']);
end

switch D
    case 5
        % not support implicit extension
        for k = 1:size(A, 3)
            A(:,:,k,:,:) = squeeze(A(:,:,k,:,:)).*mask;
        end
    otherwise
        A = A.*mask;
end

end
