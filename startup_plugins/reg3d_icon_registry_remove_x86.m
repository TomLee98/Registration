%% You can run this script if Reg3D App has been uninstalled.

% This script will remove registry items: registration related icons.

try
    % clear .regproj icon support
    [s11, ~] = system('reg delete HKEY_CLASSES_ROOT\.regproj');
    [s12, ~] = system('reg delete HKEY_CLASSES_ROOT\Reg3Dproj');

    % clear .mask icon support
    [s21, ~] = system('reg delete HKEY_CLASSES_ROOT\.mask');
    [s22, ~] = system('reg delete HKEY_CLASSES_ROOT\Reg3Dmask');

    % clear .rmv icon support
    [s31, ~] = system('reg delete HKEY_CLASSES_ROOT\.rmv');
    [s32, ~] = system('reg delete HKEY_CLASSES_ROOT\Reg3Dimage');

    if any([s11, s12, s21, s22, s31, s32])
        fprintf("Clear registry failed. Recommand check registry items:\n" + ...
            "[HKEY_CLASSES_ROOT\\.regproj]\n" + ...
            "[HKEY_CLASSES_ROOT\\Reg3Dproj]\n" + ...
            "[HKEY_CLASSES_ROOT\\.mask]\n" + ...
            "[HKEY_CLASSES_ROOT\\Reg3Dmask]\n" + ...
            "[HKEY_CLASSES_ROOT\\.rmv]\n" + ...
            "[HKEY_CLASSES_ROOT\\Reg3Dimage]\n");
    else
        fprintf("Clear registry successfully.\n");
    end
catch 
    fprintf("Clear registry failed because some unexpected error.\n")
end
