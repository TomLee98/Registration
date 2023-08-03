function status = reg3d_update(mlapp_file, auto_clean)
%REG3D_UPDATE This function auto update the reg3D app
% Input:
%   - mlapp_file: 1*1 string, the mlapp install file full path
%   - auto_clean: 1*1 bool, true for auto clean downloads folder
% Output:
%   - status: 1*1 bool, true for update script generating success, false
%             for failed

arguments
    mlapp_file (1,1) string;
    auto_clean (1,1) logical;
end

% This function generate the auto_update script 
ver_pattern = '[0-9][.][0-9][.][0-9]';
ver = regexp(mlapp_file, ver_pattern, "match");

if auto_clean == true
    clean_str = "delete(mlapp_file);";
else
    clean_str = "";
end

s = sprintf("%%This is Reg3D auto update script, which was generated automatically.\n" + ...
    "\n"+ ...
    "mlapp_file = ""%s"";\n"+...
    "matlab.apputil.install(mlapp_file);\n"+...
    "%s"+ ...
    "\ndisp(""Reg3D 更新完成！当前版本：%s"");\n", ...
    mlapp_file, clean_str, ver);

% put s in m file to userpath
file_name = fullfile(userpath, "auto_update.m");
fid = fopen(file_name,"w");
if fid == -1
    status = false;
else
    fprintf(fid,"%s",s);
    fclose(fid);
    status = true;
end

end
