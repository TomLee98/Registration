function tf_est = imregopzr(moving, fixed, rsFixed, sc_flag, tol, calg, cargs)
%IMREGOPZR This function use imregcorr and z optimization for transformation
% estimation
% Input:
%   - moving: m-by-n-by-p uint16 moving volume
%   - fixed: m-by-n-by-p uint16 reference volume
%   - rsFixed: 1-by-3 array, [x,y,z] coordinate resolution, unit as um/pix
%   - sc_flag: 1-by-1 logical, indicate if channel is structured channel,
%              true as default
%   - shift_max: 1-by-1 double, max shift at z direction, unit as pixel, 
%               2 as default
%   - tol: 1-by-1 double, z shift optimal tolerance, 1e-3 as default
%   - calg: 1-by-1 string, must be in set ["mmt","pcorr","fpp","none"],
%           "mmt" as default
%   - cargs: 1-by-1 struct, with coarse registration parameters, which
%           are ["Operator","QT","NumOctave"] at dual-channels mode and
%               ["Filter", "VT", "Radius"] as one channel mode
% Output:
%   - tf_est: 1-by-1 transltform3d object
%
% see also: imregcorr, imregmc, imregfpp, fminbnd, imwarp

arguments
    moving      (:,:,:) uint16
    fixed       (:,:,:) uint16
    rsFixed     (1,3)   double              % [x,y,z] coordinate resolution, unit as um/pix
    sc_flag     (1,1)   logical = true      % structured channel flag, true for exist
    tol         (1,1)   double {mustBeInRange(tol, 0, 1)} = 1e-3
    calg        (1,1)   string {mustBeMember(calg, ["mmt","pcorr","fpp"])} = "mmt"
    cargs       (1,1)   struct = struct("Operator","SIFT", "QT",0.0133, "NumOctave",3);
end

KUN_COEFF = 2.5;

% do maximum z projection  for imregcorr
mov_img = max(moving, [], 3);
ref_img = max(fixed, [], 3);

% get imwarp filled value
mi_val = mean([mov_img(1,:)'; mov_img(end,:)'; mov_img(:,1); mov_img(:,end)]);
ri_val = mean([ref_img(1,:)'; ref_img(end,:)'; ref_img(:,1); ref_img(:,end)]);

% add border for a square image
[height, width] = size(ref_img);
if height == width
    % already square image
    dw = [0, 0];
    dh = [0, 0];
else
    np = nextpow2(max(height, width));
    dw = [fix((2^np - width)/2), ceil((2^np - width)/2)];
    dh = [fix((2^np - height)/2), ceil((2^np - height)/2)];
end
ref_img = padarray(ref_img, [dh(1), dw(1)], "replicate","pre");
mov_img = padarray(mov_img, [dh(1), dw(1)], "replicate","pre");
ref_img = padarray(ref_img, [dh(2), dw(2)], "replicate","post");
mov_img = padarray(mov_img, [dh(2), dw(2)], "replicate","post");

% create reference coordinate
rref = imref2d(size(ref_img), rsFixed(1), rsFixed(2));
rref3d = imref3d(size(fixed), rsFixed(1), rsFixed(2), rsFixed(3));

% calculate shift on XY plane
switch calg
    case "mmt"
        % use imregmmt1 for robust shift estimation, only if one
        % predominant in the scene could be best
        if sc_flag == true
            ref_img = ref_img - ri_val;     % remove background for little objects better estimation
            mov_img = mov_img - mi_val;
            tf0 = imregmc(mov_img, ref_img, rref);
        else
            % just return the transltform3d object
            tf_est = imregmc3(moving, fixed, rref3d, cargs);
            return;
        end
    case "pcorr"
        % use imregcorr for robust shift estimation, which is
        % useful when there are high contrast scene
        % note that 'rigid' and 'affine' always return bad estimation
        tf0 = imregcorr(mov_img, rref, ref_img, rref, "translation");
    case "fpp"
        % use imregfpp for robust shift estimation, which is useful
        % when there are no strong time-varying local optical flow field
        [~, tf0] = imregfpp(mov_img, ref_img, cargs);
    otherwise
end

% estimate z0 by mess center on projection on Z axis
% as iteration initial value
Mr = reshape(sum(fixed-ri_val, [1,2]), 1, []);
Mv = reshape(sum(moving-mi_val, [1,2]), 1, []);
dz0 = sum(Mr.*(1:size(fixed,3)))/(sum(Mr)+eps) ...
    - sum(Mv.*(1:size(moving,3)))/(sum(Mv)+eps);

% optimize the z shift by immse as loss function
fminbnd_opts = optimset('MaxIter',100, 'TolX', tol);

optf = @(x)opfun(x, moving, fixed, rref3d, tf0, dz0, mi_val);

[z_, ~] = fminbnd(optf, -KUN_COEFF, KUN_COEFF, fminbnd_opts);

% omit the micro shift, which may be correction artifact
if abs(z_+dz0) < 1e-2, z_ = -dz0; end

% transform rigid2d object to affine3d object as imregtform initialized
% transformation estimation
tf_est = tformto3d(tf0, z_+dz0);

    function f = opfun(z_, mov_, ref_, ra_, tf0_, dz0_, fival_)
        T = tformto3d(tf0_, z_+dz0_);
        % imrotate3 and imtranslate?
        mov_ = imwarp(mov_, ra_, T, "linear", ...
            "FillValues",fival_, "OutputView",ra_);
        f = immse(mov_, ref_);
    end

    function T = tformto3d(tf0_, z_)
        if isa(tf0_, "affine2d")
            T = [[[tf0_.T(1:2,1:2),[0;0]];[0,0,1]];[tf0_.T(3,1:2), z_]];
            T = [T, [0;0;0;1]];
            T = affine3d(T);
        else
            T = [eye(3), [tf0_.Translation';z_]];
            T = [T; [0,0,0,1]];
            T = transltform3d(T);
        end
    end
end