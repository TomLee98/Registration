function [registered, tform] = reg_manual(template,moving,varargin)
%REGISTER2D_MAMUAL This function for manual register by using control point
% algorithm
%   Input:
%   - template: the fixed volume, 3D/4D matrix
%   - moving: the moving volume, 3D/4D matrix
%   - varargin:
%       - tform: the transformation type, could be
%               'nonreflectivesimilarity'(default),'similarity','affine',
%               'projective','polynomial','pwl','lwm'
%       - degree: parameter for pulynomial transformation, 3(default)
%       - nlwm: parameter for lwm transformation, 12(default)
%       - zproj: the z projection method, could be
%               'min','max'(default),'mean','median'
%       - loc: the location channel, could be 1(default),2,3 for R,G,B
%       - interp: the interp algorithm, could be
%                'nearest','linear','cubic'(default)
%       - smooth: the edge smooth flag, could be true, false(default)
%   Output:
%   - registered: the registered volume, 3D/4D matrix
%   - tform: the transformation estimated by fitgeotrans, which is 2D
%            transformation (on each plane)
%
%   See also fitgeotrans, imref2D, imwarp, cpselect

% REGISTER2D_MANUAL: 
% Version: 1.0.0
%   *** you can use register2D_MANUAL to align the 3D image in manual
%   *** fix the no selection bug
% 
% Version: 1.0.1
%   *** append interp option
%   *** append smooth option

% Copyright (c) 2022, Weihan Li

VALID_TRANSFORMATION = ["nonreflectivesimilarity","similarity","affine",...
    "projective","polynomial","pwl","lwm"];
VALID_ZPROJECTION = ["min","max","mean","median"];
VALID_INTERP = ["nearest","linear","cubic"];

p = inputParser;
valid_template = @(x)validateattributes(x,{'numeric'},{'nonnan'});
valid_moving = @(x)validateattributes(x,{'numeric'},{'nonnan'});
valid_tform = @(x)isscalar(x) && ismember(x,VALID_TRANSFORMATION);
valid_zproj = @(x)isscalar(x) && ismember(x,VALID_ZPROJECTION);
valid_degree = @(x)validateattributes(x,"double", ...
    {'scalar','integer','>=',2,'<=',4});
valid_nlwm = @(x)validateattributes(x,"double", ...
    {'scalar','integer','>=',6});
valid_loc = @(x)validateattributes(x,"double",{'scalar','>=',1,'<=',3});
valid_interp = @(x)isscalar(x) && ismember(x,VALID_INTERP);
valid_smooth = @(x)validateattributes(x,{'logical'},{'scalar'});

%======================== DEFAULT PARAMETER SETTING =======================
default_tform = "nonreflectivesimilarity";
default_degree = 3;
default_nlwm = 12;
default_zproj = "max";
default_loc = 1;
default_interp = "cubic";
default_smooth = false;
%==========================================================================

addRequired(p,'template',valid_template);
addRequired(p,'moving',valid_moving);
addParameter(p,'tform',default_tform,valid_tform);
addParameter(p,'degree',default_degree,valid_degree);
addParameter(p,'nlwm',default_nlwm,valid_nlwm);
addParameter(p,'zproj',default_zproj,valid_zproj);
addParameter(p,'loc',default_loc,valid_loc);
addParameter(p,'interp',default_interp,valid_interp);
addParameter(p,'smooth',default_smooth,valid_smooth);
p.parse(template,moving,varargin{:});

if ndims(p.Results.template) ~= ndims(p.Results.moving)
    error("Dimension not match!");
end

tform_type = p.Results.tform;
degree = p.Results.degree;
nlwm = p.Results.nlwm;
z_projection = p.Results.zproj;
loc = p.Results.loc;
interp = p.Results.interp;
smoothing = logical(p.Results.smooth);

registered = moving;
zip_dim = ndims(template); % because data format is [x,y,z] or [x,y,c,z]

if zip_dim == 3
    border_val = imborderval(template, 5);
elseif zip_dim == 4
    border_val = imborderval(squeeze(template(:,:,loc,:)), 5);
end
border_val = mean(border_val([1,3,5,6]));

switch z_projection
    case "min"
        template = squeeze(min(template,[],zip_dim));
        moving_zproj = squeeze(min(moving,[],zip_dim));
    case "max"
        template = squeeze(max(template,[],zip_dim));
        moving_zproj = squeeze(max(moving,[],zip_dim));
    case "mean"
        template = squeeze(mean(template,[],zip_dim));
        moving_zproj = squeeze(mean(moving,[],zip_dim));
    case "median"
        template = squeeze(median(template,[],zip_dim));
        moving_zproj = squeeze(median(moving,[],zip_dim));
    otherwise
end

% linearity scale for cpselect viewing
template = uint8(255*double(template - min(template,[],[1,2]))./...
    double(max(template,[],[1,2])-min(template,[],[1,2])));
moving_zproj = uint8(255*double(moving_zproj - min(moving_zproj,[],[1,2]))./...
    double(max(moving_zproj,[],[1,2])-min(moving_zproj,[],[1,2])));

% using cpselect to select the control points in GUI
[mp,fp] = cpselect(moving_zproj(:,:,loc),template(:,:,loc),'Wait',true);

if isempty(mp) || isempty(fp)
    tform = [];
    return;
end


switch tform_type
    case "polynomial"
        tform = fitgeotrans(mp,fp,tform_type,degree);
    case "lwm"
        tform = fitgeotrans(mp,fp,tform_type,nlwm);
    otherwise
        tform = fitgeotrans(mp,fp,tform_type);
end

% create the output view
Rfixed = imref2d(size(template(:,:,loc)));

% using imwarp apply the transform on each channel
for c = 1:size(registered,3)
    registered(:,:,c,:) = ...
            imwarp(moving(:,:,c,:),tform,interp,"SmoothEdges",smoothing,...
            "OutputView",Rfixed,'FillValues',border_val);
end

end

