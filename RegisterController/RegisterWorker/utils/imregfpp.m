function [points, tform] = imregfpp(moving, fixed, args)
% IMREGFPP: This function use feature points matching for 
% translation shift estimation
% Input:
%   - moving: m-by-n uint16 matrix, the moving image
%   - fixed: m-by-n uint16 matrix, the reference image
%   - Rf: 1-by-1 imref2d object, as reference coordinate system
% Output:
%   - points: 1-by-2 cell array, {moving, fixed}, each element is 
%     point feature object
%   - tform: 1-by-1 transltform2d object

arguments
    moving  (:,:)   uint16
    fixed   (:,:)   uint16
    args    (1,1)   struct
end

if any(size(moving)~=size(fixed))
    throw(MException("imregmmt1:imageSizeNotMatch", ...
        "Moving and fixed image must have the same size."));
end

switch args.Operator
    case "KAZE"
        ptsOriginal = detectKAZEFeatures(fixed, "Threshold",args.QT, ...
            "NumOctaves",args.NumOctave);
        ptsDistorted = detectKAZEFeatures(moving, "Threshold",args.QT, ...
            "NumOctaves",args.NumOctave);
    case "SIFT"
        ptsOriginal = detectSIFTFeatures(fixed, "ContrastThreshold",args.QT, ...
            "NumLayersInOctave",args.NumOctave);
        ptsDistorted = detectSIFTFeatures(moving, "ContrastThreshold",args.QT, ...
            "NumLayersInOctave",args.NumOctave);
    case "SURF"
        ptsOriginal = detectSURFFeatures(fixed, "MetricThreshold",args.QT, ...
            "NumOctaves",args.NumOctave);
        ptsDistorted = detectSURFFeatures(moving, "MetricThreshold",args.QT, ...
            "NumOctaves",args.NumOctave);
    case "BRISK"
        ptsOriginal = detectBRISKFeatures(fixed, "MinQuality",args.QT, ...
            "NumOctaves",args.NumOctave, "MinContrast",eps);
        ptsDistorted = detectBRISKFeatures(moving, "MinQuality",args.QT, ...
            "NumOctaves",args.NumOctave, "MinContrast",eps);
    case "Harris"
        ptsOriginal = detectHarrisFeatures(fixed, "MinQuality",args.QT);
        ptsDistorted = detectHarrisFeatures(moving, "MinQuality",args.QT);
    otherwise
        throw(MException("imregfpp:invalidDetector", ...
            "Unrecognized detector."));
end

% extract the features (with pesudo pairs)
[featuresOriginal,validPtsOriginal] = extractFeatures(fixed, ptsOriginal);
[featuresDistorted,validPtsDistorted] = extractFeatures(moving, ptsDistorted);

% match the features
index_pairs = matchFeatures(featuresOriginal,featuresDistorted,"Unique",true);

if size(index_pairs, 1) < 3
    % too less points, can not estimate the 'affine' transformation
    % return identity transformation
    points = {validPtsDistorted, validPtsOriginal};
    tform = transltform2d();
    return;
end

% extract the matched point pairs
matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));

% estimate the 2d translation transformation
try
    [tform,inlierIdx] = ...
        estgeotform2d(matchedPtsDistorted, matchedPtsOriginal, ...
        "affine", "Confidence",99, "MaxDistance",2.5);

    tform = tfcast(tform);

    inlierPtsOriginal  = matchedPtsOriginal(inlierIdx,:);
    inlierPtsDistorted = matchedPtsDistorted(inlierIdx,:);

    points = {inlierPtsDistorted, inlierPtsOriginal};
catch
    % estimation with strong constrains so that there is matching conflict
    % between matchFeatures and estimateGeometricTransform2D/estgeotform2d
    warning("imregfpp:fewerMatchingPoints", ...
        "No enough points matching, coarse estimation may be bad.");
    points = {validPtsDistorted, validPtsOriginal};
    tform = transltform2d();
end

end

% ======================= utilities function ===========================

function A = tfcast(A)
% This function cast the transformation to translation only
if isa(A, "affinetform2d")
    % remove shear and rotation
    A.A(1:2, 1:2) = eye(2);
    A = transltform2d(A.A);
else
    A = transltform2d();
end
end
