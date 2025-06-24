function [tf_est, movol_est] = imregcoarse(moving, fixed, rsFixed, sc_flag, tol, calg, cargs)
%IMREGCOARSE This function use multiple robust coarse registration
%algorithm for violently motion correction
% Input:
%   - moving: m-by-n-by-p uint16 moving volume
%   - fixed: m-by-n-by-p uint16 reference volume
%   - rsFixed: 1-by-3 array, [x,y,z] coordinate resolution, unit as um/pix
%   - sc_flag: 1-by-1 logical, indicate if channel is structured channel,
%              true as default
%   - tol: 1-by-1 double, z shift optimal tolerance, 1e-3 as default
%   - calg: 1-by-1 string, must be in set ["mmt","pcorr","fpp","none"],
%           "mmt" as default
%   - cargs: 1-by-1 struct, with coarse registration parameters, which
%           are ["Operator","QT","NumOctave"] at dual-channels mode and
%               ["Filter", "VT", "Radius"] as one channel mode
% Output:
%   - tf_est: 1-by-1 transltform3d object
%   - movol_est: m-by-n-by-p uint16 registered array
%
% see also: imregcorr, imregmc, imregfpp, fminbnd, imwarp
arguments
    moving      (:,:,:) uint16
    fixed       (:,:,:) uint16
    rsFixed     (1,3)   double      % [x,y,z] coordinate resolution, unit as um/pix
    sc_flag     (1,1)   logical = true
    tol         (1,1)   double {mustBeInRange(tol, 0, 1)} = 1e-3
    calg        (1,1)   string {mustBeMember(calg, ["mmt","pcorr","fpp","none"])} = "mmt"
    cargs       (1,1)   struct {mustBeCoarseRegistrationArguments} = ...
                        struct("Operator","SIFT", "QT",0.0133, "NumOctave",3);
end

switch calg
    case "none"
        % only 3d identity transformation
        tf_est = transltform3d();
    otherwise
        % pseudo 3d transformation: XY translation + Z translation
        tf_est = imregopzr(moving, fixed, rsFixed, sc_flag, tol, calg, cargs);
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

if ~isempty(setxor(fields_, ["Operator", "QT", "NumOctave"])) && ...
        ~isempty(setxor(fields_, ["Filter","VT","Radius"]))
    throw(MException("mustBeCoarseRegistrationArguments:invalidArgumentsItem", ...
        "Unrecognized aeguments."));
end
end
