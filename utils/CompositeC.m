function C = CompositeC(A, B, method, order, pseudoColor)
% COMPOSITEC This function composites the color channel and genetates the RGB
% representaton
% Input:
%   - A: 2D/3D uint8 matrix for channel A
%   - B: 2D/3D uint8 matrix for channel B
%   - method: the combination method, which could be "none", "red-cyan" or
%   "green-magenta", where the format follows channel "A-B"
%   - order: if no combiantion, order deteremine which channel is selected,
%   and 1 for channel A, 2 for channel B
%   - pseudoColor: the pseudo color need to viewing, which could be "r",
%   "g","mix", where "r" for "red"(~694nm light wave), and "g" for
%   "green"(~521nm light wave), "mix" for dual color viewing
% Output:
%   - C: 3D/4D uint8 true color matrix

% Version 1.1.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A (:,:,:) {isnumeric(A)}
    B (:,:,:) {isnumeric(B)}
    method (1,1) string {ismember(method,["none","red-cyan",...
        "green-magenta"])} = "none";
    order (1,1) double {ismember(order,[1,2])} = 1;
    pseudoColor (1,1) string {ismember(pseudoColor, ["r","g","mix"])} = "r";
end

if ~isempty(A)&&~isempty(B)&&any(size(A)~=size(B))
    error("Size not match!");
end

switch method
    case "none"
        % no combination, only single pseudo color
        if order == 1, C = A; else, C = B; end
        switch pseudoColor
            case "r"
                % generate the pure red for 694nm light wave
                % on both image and viewer3d solution
                C = cat(ndims(C)+1, C, zeros([size(C),2],'like',C));
            case "g"
                if ismatrix(C)
                    % generate the mixing green for 521nm light wave
                    % image view solution need both R and G
                    C = cat(ndims(C)+1, uint8(0.4*double(C)),C,zeros(size(C),'like',C));
                else
                    % generate the mixing green for 521nm light wave
                    % viewer3d only need G channel
                    C = cat(ndims(C)+1, zeros(size(C),'like',C),C,zeros(size(C),'like',C));
                end
            case "mix"
                % generate the red and green mixing view
                % B for red channel and A for green, where using mixer with
                % (0.8, 0.2) for color tuning
                C = cat(ndims(C)+1,uint8(0.8*double(B)+0.08*double(A)),A,zeros(size(C),'like',C));
            otherwise
                % gray scale(RGB representation)
                % comments: we will not going here
                C = cat(ndims(C),C,C,C);
        end
    case "red-cyan"
        C = cat(ndims(A)+1,A,B,B);
    case "green-magenta"
        C = cat(ndims(A)+1,B,A,B);
end
end

