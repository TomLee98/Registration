function tf_est = imregopzr(moving, fixed, rsFixed, shift_max, tol, calg, cargs)
%IMREGOPZR This function use imregcorr and z optimization for transformation
% estimation on platform version < R2022b
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
%   - tf_est:
%
% see also: imregcorr, imregmc, imregfpp, fminbnd, imwarp

arguments
    moving      (:,:,:) uint16
    fixed       (:,:,:) uint16
    rsFixed     (1,3)   double      % [x,y,z] coordinate resolution, unit as um/pix
    shift_max   (1,1)   double {mustBeNonnegative} = 2
    tol         (1,1)   double {mustBeInRange(tol, 0, 1)} = 1e-3
    calg        (1,1)   string {mustBeMember(calg, ["mmt","pcorr","fpp"])} = "mmt"
    cargs       (1,1)   struct = struct("Operator","SIFT", "QT",0.0133, "NumOctave",3);
end

zlim = [-1, 1]*shift_max;

% do maximum z projection  for imregcorr
mov_img = max(moving, [], 3);
ref_img = max(fixed, [], 3);

% get imwarp filled value
fi_val = mean([mov_img(1,:)'; mov_img(end,:)'; mov_img(:,1); mov_img(:,end)]);

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
        tf0 = imregmc(mov_img, ref_img, rref);
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

% optimize the z shift by immse as loss function
fminbnd_opts = optimset('MaxIter',100, 'TolX', tol);

optf = @(x)opfun(x, moving, fixed, rref3d, tf0, fi_val);

[z_, ~] = fminbnd(optf, zlim(1), zlim(2), fminbnd_opts);

% omit the micro shift, which may be correction artifact
if abs(z_) < 1e-2, z_ = 0; end

% transform rigid2d object to affine3d object as imregtform initialized
% transformation estimation
tf_est = tformto3d(tf0, z_);

    function f = opfun(z_, mov_, ref_, ra_, tf0_, fival_)
        T = tformto3d(tf0_, z_);
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