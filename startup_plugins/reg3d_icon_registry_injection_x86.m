function status = reg3d_icon_registry_injection_x86(temp_folder)
%REG3D_ICON_INJECTION This function injects icon which connects with
%*.regproj file on PC platform
arguments
    temp_folder (1,1)   string  {mustBeFolder}
end

% get icon file full path
rf_icon = mfilename("fullpath");
[rf_icon, ~, ~] = fileparts(rf_icon);
rf_icon = string(rf_icon).replace("startup_plugins", "sources") + filesep + "icon.png";
rf_icon = rf_icon.replace("\", "\\");

%% create .reg file (register file on Windows: Version 5.00)
%[TEXT START]
rfp_file = fullfile(temp_folder, "regicon.reg");

% clear temporary file
status = remove_temporary_file(rfp_file);
if status==-1, return; end

try
    fid = fopen(rfp_file, "a+", "n", "UTF-8");

    %Windows Registry Editor Version 5.00
    fprintf(fid, "Windows Registry Editor Version 5.00\n\n");

    %[HKEY_CLASSES_ROOT\.regproj]
    %@="Reg3Dfile"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\.regproj]\n");
    fprintf(fid, "@=""Reg3Dfile""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dfile]
    %@="Registration Project File"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dfile]\n");
    fprintf(fid, "@=""Registration Project File""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dfile\DefaultIcon]
    %@="<folder>\*.png"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dfile\\DefaultIcon]\n");
    fprintf(fid, "@=""%s""\n\n", rf_icon);

    %[HKEY_CLASSES_ROOT\Reg3Dfile\shell]
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dfile\\shell]\n");

    %[TEXT END]
    fclose(fid);
catch
    warning("registry_injection:invalidPermission", ...
        "Temporary registry file creation was denied.");
    status = -1;
    return;
end

%% inject by system call
cmd = ['reg import ', char(rfp_file)];

try
    % This command need to run as "Administrator"
    [status, ~] = system(cmd);
    if status ~= 0
        warning("registry_injection:invalidPermission", ...
            "Registry changes imported failed. You could import manually as administrator.");
        fprintf("Regustry file location: %s\n", rfp_file);
    else
        % clear temporary file
        status = remove_temporary_file(rfp_file);

        fprintf("Registry changes imported successfully.\n" + ...
            "They will take effect after explorer.exe restart.\n");
    end
catch ME
    warning("registry_injection:unknownException", ...
        "An exception was thrown: %s.\n", ME.identifier);
    status = -1;
    return;
end

end

function status = remove_temporary_file(rfp_file)
try
    delete(rfp_file);   % remove the old file
    status = 0;
catch
    warning("registry_injection:invalidPermission", ...
        "Cleaning temporary registry file was denied.");
    status = -1;
    return;
end
end