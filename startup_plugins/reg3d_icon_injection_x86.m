function status = reg3d_icon_injection_x86(temp_folder)
%REG3D_ICON_INJECTION This function injects icon which connects with
%*.regproj file on PC platform
arguments
    temp_folder (1,1)   string  {mustBeFolder}
end

% get icon file full path
rf_icon = mfilename("fullpath");
[rf_icon, ~, ~] = fileparts(rf_icon);
rf_icon = string(rf_icon).replace("startup_plugins", "sources") + filesep + "icon.png";

%% create .reg file (register file on Windows: Version 5.00)
%[TEXT START]
rfp_file = fullfile(temp_folder, "regicon.reg");
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

%% inject by system call

cmd = ['regedit.exe /s "', rfp_file, '"'];

try
    [status, ~] = system(cmd);
    if status ~= 0
        fprintf("Registry changes imported failed.\n");
    else
        fprintf("Registry changes imported successfully. They will take effect after a restart.\n");
    end
catch ME
    fprintf("An exception was thrown: %s", ME.identifier);
    status = -1;
end

end

