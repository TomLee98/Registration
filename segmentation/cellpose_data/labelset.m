% This script reads all slices from ./trainset and generates mask by
% manual labeling
% Note: 
%%
clear;
clc;

%% set parameters
auto_contrast = true;

%%
pathToThis = fileparts(mfilename('fullpath'));
src_folder = [pathToThis, '\trainset'];

list = string({dir(src_folder).name});

ispng = list.endsWith(".png","IgnoreCase",true);

if any(ispng)
    img_files = string(src_folder) + "\" + list(ispng);
    img_files = img_files(img_files.endsWith("_im.png") ...
        & ~isfile(img_files.replace("_im.png", "_mask.png")));
    
    % Man in loop
    for k = 1:numel(img_files)
        %% check if mask file is exist
        if ~img_files(k).endsWith("_im.png"), continue, end
        mask_file = img_files(k).extractBefore("_im.png") + "_mask.png";
        if exist(mask_file, "file"), continue; end

        [~, name, ext] = fileparts(img_files(k));
        fig = figure("Name","<Space> for next ROI, <Enter> for submit. ");
        
        set(fig, 'windowkeypressfcn', @keypressfcn);
        fig.UserData = mask_file;
        
        img = imread(img_files(k));
        h = imshow(uint8(rescale(img, 0, 255)));

        setFigure(fig, [0,0,0.5,0.8]);

        title(sprintf("Name:%s     Progress:%d/%d", name+ext, k, numel(img_files)), ...
            "FontSize",16, "Interpreter","none");

        waitfor(h);
    end

    fprintf("All labels(%d) well done.", numel(img_files));
end

%%
function keypressfcn(h, evt)
if evt.Character == sprintf("\r")
    % export ROI and close handle
    export_mask(h);
    
    delete(h);
elseif evt.Character == " "
    % draw new ROI
    drawfreehand();
end
end

% This function exports mask from ROIs
function export_mask(h)
export_file = h.UserData;
ax = findobj(h, "Type", "Axes");
mask = false(floor(ax.YLim(end)), floor(ax.XLim(end)));

rois = findobj(ax.Children, "Type", "images.roi.Freehand");

for n = 1:numel(rois)
    mask = mask | createMask(rois(n));
end

imwrite(mask, export_file);

end