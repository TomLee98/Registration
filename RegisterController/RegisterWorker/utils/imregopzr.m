function [tf_est, movol_est] = imregopzr(moving, fixed, rsFixed, shift_max, tol, calg)
%IMREGOPZR This function use imregcorr and z optimization for transformation
% estimation on platform version < R2022b
% Input:
%   - moving:
%   - fixed:
%   - rsFixed:
%   - tformType:
%   - shift_max:
%   - tol:
%   - calg: 
% Output:
%   - tf_est:
%   - movol_est:
%
% see also: imregcorr, fminbnd, imwarp
arguments
    moving      (:,:,:) uint16
    fixed       (:,:,:) uint16
    rsFixed     (1,3)   double      % [x,y,z] coordinate resolution, unit as um/pix
    shift_max   (1,1)   double {mustBeNonnegative} = 2
    tol         (1,1)   double {mustBeInRange(tol, 0, 1)} = 1e-3
    calg        (1,1)   string {mustBeMember(calg, ["mmt","pcorr","both","none"])} = "mmt"
end

% compatibility flag
is_matlab_old_ver = isMATLABReleaseOlderThan("R2022b");

if calg == "none"
    if is_matlab_old_ver == true
        tf_est = affine3d();
    else
        tf_est = transltform3d();
    end
    movol_est = moving;
    return;
end

% zlim = size(fixed, 3)*[-1, 1];
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

% ================= (optional) first step branch =====================
if (calg == "mmt") || (calg == "both")
    % use imregmmt1 for the most robust shift estimation
    tf0 = imregmc(mov_img, ref_img, rref);

    % use tf0 do first step registration
    mov_img = imwarp(mov_img, rref, tf0, "linear", "FillValues", fi_val);
else
    % identity transformation instead
    if is_matlab_old_ver == true
        tf0 = affine2d();
    else
        tf0 = transltform2d();
    end
end

% ================= (optional) second step branch =====================
if (calg == "pcorr") || (calg == "both")
    % 2-D rigid transformation estimation
    % tf2d_: affine2d or transltform2d
    tf2d_ = imregcorr(mov_img, rref, ref_img, rref, "translation");  
else
    % identity transformation instead
    if is_matlab_old_ver == true
        tf2d_ = affine2d();
    else
        tf2d_ = transltform2d();
    end
end


% ================ (acquired) z shift optimal =======================
fminbnd_opts = optimset('MaxIter',100, 'TolX', tol);

% optimize the z shift by immse as loss function
optf = @(x)opfun(x, moving, fixed, rref3d, tf0, tf2d_, fi_val);
[z_, ~] = fminbnd(optf, zlim(1), zlim(2), fminbnd_opts);

% omit the micro shift, which may be correction artifact
if abs(z_) < 1e-2, z_ = 0; end

% transform rigid2d object to affine3d object as imregtform initialized
% transformation estimation
tf_est = tformto3d(tf0, tf2d_, z_);

if nargout == 2
    % 1 memory copy from <imwarp>
    movol_est = imwarp(moving, rref3d, tf_est, "linear",...
        "OutputView",rref3d, 'FillValues',fi_val);
end

    function f = opfun(z_, mov_, ref_, ra_, tf0_, tf_, fival_)
        T = tformto3d(tf0_, tf_, z_);
        % imrotate3 and imtranslate?
        mov_ = imwarp(mov_, ra_, T, "linear", ...
            "FillValues",fival_, "OutputView",ra_);
        f = immse(mov_, ref_);
    end

    function T = tformto3d(tf0_, tf_, z_)
        if isa(tf_, "affine2d")
            T = [[[tf_.T(1:2,1:2),[0;0]];[0,0,1]];[tf_.T(3,1:2)+tf0_.T(3,1:2), z_]];
            T = [T, [0;0;0;1]];
            T = affine3d(T);
        else
            T = [eye(3), [(tf_.Translation+tf0_.Translation)';z_]];
            T = [T; [0,0,0,1]];
            T = transltform3d(T);
        end
    end
end