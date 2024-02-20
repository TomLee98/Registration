function [status, url, ver_latest] = acquire_update(url, ver)
%ACQUIRE_UPDATE This function compare localhost app version and cloud version,
% return the update status and source url
% Input:
%   - ver: 1*1 string, the localhost app version
% Oupput:
%   - status: 1*1 string, the updatable flag. "UPDATABLE" for updatable, 
%           "NO_UPDATE" for no new version, "ACCESS_FAILED" for remote 
%           repository accessing failed
%   - url: char array, the new app source url

% Version: 1.0.0
% VERSION DEFINATION: <MAIN_VERSION>.<SUB_VERSION>.<PATCH_VERSION>

% Copyright (c) 2023, Weihan Li

arguments
    url (1,1) string;
    ver (1,1) string;
end

% read the repository and check the version
ver_pattern = '[0-9][.][0-9][.][0-9]';

listing = dir(url);
if isempty(listing)
    % NAS file server
    if lower(string(url)).contains("nas")
        status = "ACCESS_FAILED";
        url = "";
        return;
    end
end

files = struct2cell(listing);
files_mlapp = string(files(1,:));

ver_repos =  regexp(files_mlapp,ver_pattern,"match");
ver_repos(cellfun(@isempty, ver_repos)) = [];

% transform version to number for comparing
ver_repos =  fliplr(string(ver_repos));
ver_latest = latest_version(ver_repos);

if compare_version(ver, ver_latest) == -1
    status = "UPDATABLE";
    url_idx = files_mlapp.contains(ver_latest);
    url = files(1:2, url_idx);
    url = fullfile(url{2}, url{1});
else
    status = "NO_UPDATE";
    url = "";
end

ver_latest = ver_latest.char();
end

function res = compare_version(ver_lhs, ver_rhs)
% this function return the version status between ver1 and ver2
% Input:
%   - ver_lhs: string, the left hand side operation version
%   - ver_rhs: string, the right hand side operation version
% Output:
%   - res: logical, the left version level relative to right version

ver_lhs_arr = str2double(ver_lhs.split("."));
ver_rhs_arr = str2double(ver_rhs.split("."));

if all(ver_lhs_arr==ver_rhs_arr)
    res = 0;
else
    if ver_lhs_arr(1)<ver_rhs_arr(1)
        res = -1;
    elseif ver_lhs_arr(1)>ver_rhs_arr(1)
        res = 1;
    else
        if ver_lhs_arr(2)<ver_rhs_arr(2)
            res = -1;
        elseif ver_lhs_arr(2)>ver_rhs_arr(2)
            res = 1;
        else
            if ver_lhs_arr(3)<ver_rhs_arr(3)
                res = -1;
            else
                res = 1;
            end
        end
    end
end

end

function ver_latest = latest_version(vers)
%   this function get the latest version in vers
% Input:
%   - vers: string vector
% Output:
%   - ver_latest: 1*1 string, the latest version in vers

arguments
    vers string {isvector};
end

ver_latest = "1.0.0";
for k = 1:numel(vers)
    if compare_version(ver_latest, vers(k)) == -1
        ver_latest = vers(k);
    end
end
end
