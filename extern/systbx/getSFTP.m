function sftp_handle = getSFTP(host)
%GETSFTP This function call LabNasLogin.mlapp and will be protected by
%pcode for LabNas user name and password security
% Assume that pcode inner is protected

sftp_handle = [];
N = 3;

lnl_handle = LabNasLogin();
waitfor(lnl_handle,"login_status",true);

if ~isvalid(lnl_handle)
    return;
end

% take the information and delete app
username = lnl_handle.username;
if ispc()
    aes_unlocker = aesobj();
    password = aes_unlocker.decrypt(lnl_handle.password);
elseif isunix()
    % security warning:
    % may lost password, bugs on unix for AES object
    password = lnl_handle.password;
end

% try to connect the labnas with sftp
n = 0;
while isempty(sftp_handle) && (n < N)
    try
        sftp_handle = sftp(host, username, ...
            "Password",password,"ServerSystem","unix");
    catch exception
        switch exception.identifier
            case 'MATLAB:io:ftp:ftp:BadLogin'
                warning("LabNasLogin:invalidInputs", ...
                    "Invalid user name or password.");
            case 'MATLAB:io:ftp:ftp:NoConnection'
                warning("LabNasLogin:remoteNoResponse", ...
                    "[%d]Remote host no response, try to reconnect...",n+1);
            case 'MATLAB:io:ftp:ftp:InternalError'
                warning("LabNasLogin:ftpInternalError", ...
                    "[%d]Remote host no response, try to reconnect...",n+1);
            otherwise
                throwAsCaller(exception);
        end
    end
    n = n + 1;
end

delete(lnl_handle);

end

