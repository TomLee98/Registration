function A = Projection(A, method, dim, range)
% PROJECTION This function decrease dimension at the target dimension by
% using the method on array A
% Input:
%   - A: the gray image stack
%   - method: the projection method, which could be "max", "min",
%     "median", "mean"
%   - dim: the dimension need to be projection
%   - range: the data selected range need to be processed
% Output:
%   - A: the projected array 

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A (:,:,:,:,:);
    method (1,1) {ismember(method,["max","min","median","mean"])} = "mean";
    dim (1,1) double {mustBeInRange(dim,1,5)} = 4;
    range {mustBeTextScalar} = string(1:size(A, dim)).join(",");
end

% calculate the slice of A
eval_str = "(:,:,:,:)";
if dim >= 1 && dim <= 4
    insert_str = "["+range+"]" + ",";
    eval_str = insertAfter(eval_str, 2*dim-1, insert_str);
    A = eval("A"+eval_str);
else
    A = A(:,:,:,:,str2num(range)); %#ok<ST2NM> 
end

method_func = str2func(method);
if ismember(method,["mean","median"])
    A = squeeze(method_func(A,dim));
elseif ismember(method,["max","min"])
    A = squeeze(method_func(A,[],dim));
end

% recover the color channel
if ndims(A) == 3 && dim > 3
    A = reshape(A,[size(A,[1,2]),1,size(A,3)]);
end

end