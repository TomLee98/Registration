function fixvol = GetFixedVol(movie, fixed, bigflag)
% GETFIXEDVOL This function extract the fixed volume by using fixed
% information
% Input:
%   - movie: the image data, which could be high dimension array in memory
%     or a file pointer(matfile type)
%   - fixed: 1*2 string cell with {method, range}, where method could be
%     "mean", "median", "min", "max", and range is a string with positive 
%     integer numbers as its information
%   - bigflag: the big file indicator, which could be "on" or "off", "off"
%     as default
% Output:
%   - fixvol: the final fixed volume (dimension: xyz)

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    movie;
    fixed (1,2) string;   % fixed is cell, with {method, range}
    bigflag (1,1) string {ismember(bigflag,["on","off"])} = "off";
end

range = str2num(fixed(2)); %#ok<ST2NM>

if bigflag == "on"
    varlist = who(movie);
    size_mv = size(movie,varlist{1});
    fixvol = zeros(size_mv(1:end-1),'like',movie.(varlist{1})(:,:,:,1));
    switch fixed(1)
        case 'mean'
            for page = range
                fixvol = fixvol+movie.(varlist{1})(:,:,:,page);
            end
            fixvol = fixvol/numel(range); %#ok<*NASGU> 
        case 'max'
            for page = range
                fixvol = max(fixvol, movie.(varlist{1})(:,:,:,page));
            end
        case 'min'
            fixvol = inf(size_mv(1:end-1));
            for page = range
                fixvol = min(fixvol,movie.(varlist{1})(:,:,:,page));
            end
        case 'median'
            fixvol = median(movie.(varlist{1}), 4);
        otherwise
    end
else
    switch fixed(1)
        case 'mean'
            fixvol = mean(movie(:,:,:,range),4);
        case 'median'
            fixvol = median(movie(:,:,:,range),4);
        case 'max'
            fixvol = max(movie(:,:,:,range),[],4);
        case 'min'
            fixvol = min(movie(:,:,:,range),[],4);
        otherwise
    end
end

% retype fixvol from double to type of movie
switch bigflag
    case "on"
        fixvol = eval([class(movie.(varlist{1})(:,:,:,1)),'(fixvol)']);
    case "off"
        fixvol = eval([class(movie),'(fixvol)']);
end

end

