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
%   - tform: 1-by-1 affine2d(MATLAB<R2022b) or transltform2d object

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

    case "SURF"

    case "BRISK"

    case "Harris"

    otherwise
end

[featuresOriginal,validPtsOriginal] = extractFeatures(original,ptsOriginal);
[featuresDistorted,validPtsDistorted] = extractFeatures(distorted,ptsDistorted);

index_pairs = matchFeatures(featuresOriginal,featuresDistorted);

matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));

if isMATLABReleaseOlderThan("R2022b")
    [tform,inlierIdx] = estimateGeometricTransform2D(matchedPtsDistorted, matchedPtsOriginal, "affine");
else
    [tform,inlierIdx] = estgeotform2d(matchedPtsDistorted, matchedPtsOriginal, "affine");
end

inlierPtsOriginal  = matchedPtsOriginal(inlierIdx,:);
inlierPtsDistorted = matchedPtsDistorted(inlierIdx,:);

points = {inlierPtsDistorted, inlierPtsOriginal};
end

