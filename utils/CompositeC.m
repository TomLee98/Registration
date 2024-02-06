function C = CompositeC(A, B, method, order, pseudoColor)
% COMPOSITEC This function composites the color channel and genetates the RGB
% representaton
% Input:
%   - A: 2D uint8 matrix for channel A
%   - B: 2D uint8 matrix for channel B
%   - method: the combination method, which could be "none", "red-cyan" or
%   "green-magenta", where the format follows channel "A-B"
%   - order: if no combiantion, order deteremine which channel is selected,
%   and 1 for channel A, 2 for channel B
%   - pseudoColor: the pseudo color need to viewing, which could be "r",
%   "g","mix", where "r" for "red"(~694nm light wave), and "g" for
%   "green"(~521nm light wave), "mix" for dual color viewing
% Output:
%   - C: 3D uint8 true color matrix

% Version 1.1.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A (:,:) {isnumeric(A)}
    B (:,:) {isnumeric(B)}
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
                C = cat(3,C,zeros([size(C),2],'like',C));
            case "g"
                % generate the mixing green for 521nm light wave
                C = cat(3,uint8(0.4*double(C)),C,zeros(size(C),'like',C));
            case "mix"
                % generate the red and green mixing view
                % B for red channel and A for green, where using mixer with
                % (0.8, 0.2) for color tuning
                C = cat(3,uint8(0.8*double(B)+0.08*double(A)),A,zeros(size(C),'like',C));
            otherwise
                % gray scale(RGB representation)
                % comments: we will not going here
                C = cat(3,C,C,C);
        end
    case "red-cyan"
        C = cat(3,A,B,B);
    case "green-magenta"
        C = cat(3,B,A,B);
end
end

