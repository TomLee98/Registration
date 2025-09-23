function status = reg3d_icon_registry_injection_x86(temp_folder)
%REG3D_ICON_INJECTION This function injects icon which connects with
%*.regproj file on PC platform
arguments
    temp_folder (1,1)   string  {mustBeFolder}
end

% get icon file full path
[mfolder, ~, ~] = fileparts(mfilename("fullpath"));
icons_folder = string(mfolder).replace("startup_plugins", "sources");
rf_proj_icon = replace(icons_folder + filesep + constdef.PROJECT_FILE_ICON, "\", "\\");
rf_mask_icon = replace(icons_folder + filesep + constdef.MASK_FILE_ICON, "\", "\\");
rf_rmv_icon = replace(icons_folder + filesep + constdef.REGMOV_FILE_ICON, "\", "\\");

%% create .reg file (register file on Windows: Version 5.00)
%[TEXT START]
rfp_file = fullfile(temp_folder, "regicons.reg");

% clear temporary file
status = remove_temporary_file(rfp_file);
if status==-1, return; end

try
    fid = fopen(rfp_file, "a+", "n", "UTF-8");

    %% TITLE
    %Windows Registry Editor Version 5.00
    fprintf(fid, "Windows Registry Editor Version 5.00\n\n");

    %% SUPPORT FOR *.REGPROJ FILE (PROJECT FILE)
    %[HKEY_CLASSES_ROOT\.regproj]
    %@="Reg3Dproject"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\.regproj]\n");
    fprintf(fid, "@=""Reg3Dproject""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dproject]
    %@="Registration Project File"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dproject]\n");
    fprintf(fid, "@=""Registration Project File""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dproject\DefaultIcon]
    %@="<folder>\*.png"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dproject\\DefaultIcon]\n");
    fprintf(fid, "@=""%s""\n\n", rf_proj_icon);

    %[HKEY_CLASSES_ROOT\Reg3Dproject\shell]
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dproject\\shell]\n\n");

    %% SUPPORT FOR *.MASK FILE (SEGMENTOR->MASK FILE)
    %[HKEY_CLASSES_ROOT\.mask]
    %@="Reg3Dfile"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\.mask]\n");
    fprintf(fid, "@=""Reg3Dmask""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dfile]
    %@="Registration Project File"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dmask]\n");
    fprintf(fid, "@=""Registration Mask File""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dmask\DefaultIcon]
    %@="<folder>\*.png"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dmask\\DefaultIcon]\n");
    fprintf(fid, "@=""%s""\n\n", rf_mask_icon);

    %[HKEY_CLASSES_ROOT\Reg3Dfile\shell]
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dmask\\shell]\n\n");

    %% SUPPORT FOR *.RMV FILE (REGMOV FILE)
    %[HKEY_CLASSES_ROOT\.rmv]
    %@="Reg3Dimage"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\.rmv]\n");
    fprintf(fid, "@=""Reg3Dimage""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dimage]
    %@="Registration Image File"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dimage]\n");
    fprintf(fid, "@=""Registration Image File""\n\n");

    %[HKEY_CLASSES_ROOT\Reg3Dimage\DefaultIcon]
    %@="<folder>\*.png"
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dimage\\DefaultIcon]\n");
    fprintf(fid, "@=""%s""\n\n", rf_rmv_icon);

    %[HKEY_CLASSES_ROOT\Reg3Dimage\shell]
    fprintf(fid, "[HKEY_CLASSES_ROOT\\Reg3Dimage\\shell]\n");

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