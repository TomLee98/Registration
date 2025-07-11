function tf = isPyReady()
%ISPYREADY This function returns python environment status
% check the environment

% Versions of Python Compatible with MATLAB Products by Release
% https://ww2.mathworks.cn/support/requirements/python-compatibility.html
MAT_PY_VALID_MAPPING = struct("R2022b", ["2.7","3.8","3.9","3.10"], ... ...
                              "R2023a", ["3.8","3.9","3.10"], ...
                              "R2023b", ["3.9","3.10","3.11"], ...
                              "R2024a", ["3.9","3.10","3.11"], ...
                              "R2024b", ["3.9","3.10","3.11","3.12"], ...
                              "R2025a", ["3.9","3.10","3.11","3.12"]);

if isMATLABReleaseOlderThan("R2022b")
    tf = false;
    warning("isPyReady:tooOldMATLABVersion", "Fase tiff IO based on python " + ...
        "need your MATLAB version >= R2022b");
    return;
end

% get matlab version
mat_ver = string(version).extractBetween("(",")");

pe = pyenv;

if (pe.Version ~= "")
    if ismember(pe.Version, MAT_PY_VALID_MAPPING.(mat_ver))
        if ispc()
            % redirect output to avoid console display
            s_np = system("pip list | findstr numpy > ~.sftxt");
            s_tf  = system("pip list | findstr tifffile >> ~.sftxt");
        elseif isunix()
            s_np = system("pip list | grep numpy > ~.sftxt");
            s_tf  = system("pip list | grep tifffile >> ~.sftxt");
        else
            throw(MException("isPyReady:invalidOperationSystem", ...
                "Sorry, your operation system has not been supported."));
        end
        delete("~.sftxt");    % remove the redirect temporary file
        if (s_np==0) && (s_tf==0)
            % there are python and tifffile package, use tifffile replace
            % bio-formats java package
            tf = true;
        else
            tf = false;
            warning("savetiff:noTFOrNP", ...
                "Fast saving disabled because of numpy or tifffile package(s) lost.");
        end
    else
        tf = false;
        warning("savetiff:unmatchedPythonVersion","Fast saving disabled because of " + ...
            "python version is not match.");
    end
else
    tf = false;
    warning("savetiff:noPythonOrNotPC","Fast saving disabled because of " + ...
        "unsupported platform or python environment is not installed.");
end

end

