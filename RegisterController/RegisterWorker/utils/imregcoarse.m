function [tf_est, movol_est] = imregcoarse(moving, fixed, rsFixed, tfType, shift_max, tol, calg, cargs)
%IMREGCOARSE This function use multiple robust coarse registration
%algorithm for violently motion correction
% Input:
%   - moving:
%   - fixed:
%   - rsFixed:
%   - tfType: 
%   - shift_max:
%   - tol:
%   - calg: 
%   - cargs:
% Output:
%   - tf_est: 1-by-1 transltform2d, rigidtform2d or affinetform2d object
%   - movol_est:
%
% see also: imregcorr, imregmc, imregfpp, fminbnd, imwarp
arguments
    moving      (:,:,:) uint16
    fixed       (:,:,:) uint16
    rsFixed     (1,3)   double      % [x,y,z] coordinate resolution, unit as um/pix
    tfType      (1,1)   string {mustBeMember(tfType, ["translation","rigid","affine"])} = "translation"
    shift_max   (1,1)   double {mustBeNonnegative} = 2
    tol         (1,1)   double {mustBeInRange(tol, 0, 1)} = 1e-3
    calg        (1,1)   string {mustBeMember(calg, ["mmt","pcorr","fpp","uepp","none"])} = "mmt"
    cargs       (1,1)   struct {mustBeCoarseRegistrationArguments} = ...
                        struct("Operator","SIFT", "QT",0.0133, "NumOctave",3);
end

switch calg
    case "none"
        % only 3d identity transformation
        tf_est = transltform3d();
    case "uepp"
        % real 3d transformation
        tf_est = imreguepp(moving, fixed, rsFixed, tfType);
    otherwise
        % pseudo 3d transformation: XY translation + Z translation
        tf_est = imregopzr(moving, fixed, rsFixed, shift_max, ...
            tol, calg, cargs);
end

if nargout == 2
    rref3d = imref3d(size(fixed), rsFixed(1), rsFixed(2), rsFixed(3));
    % do maximum z projection  for imregcorr
    mov_img = max(moving, [], 3);

    % get imwarp filled value
    fi_val = mean([mov_img(1,:)'; mov_img(end,:)'; mov_img(:,1); mov_img(:,end)]);

    % 1 memory copy from <imwarp>
    movol_est = imwarp(moving, rref3d, tf_est, "linear",...
        "OutputView",rref3d, 'FillValues',fi_val);
end

end

function mustBeCoarseRegistrationArguments(A)
fields_ = fieldnames(A);

if ~isempty(setxor(fields_, ["Operator", "QT", "NumOctave"]))
    throw(MException("mustBeCoarseRegistrationArguments:invalidArgumentsItem", ...
        "Unrecognized aeguments."));
end
end
