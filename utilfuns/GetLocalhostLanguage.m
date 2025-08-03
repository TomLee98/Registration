function lang = GetLocalhostLanguage()
%GETREGIONINFO This function get the region setting by using systeminfo
% at Windows platform or locale at unix platform
% Input:
% Output:
%   - lang: 1-by-1 string, with localhost language setting

if ispc()
    [~, sysinfo] = system("systeminfo");
    sysinfo = lower(string(sysinfo));
    if sysinfo.contains("zh-cn")
        lang = "zh_CN";
    else
        lang = "en_US";
    end
else
    [~, sysinfo] = system("locale");
    sysinfo = lower(string(sysinfo).split("="));
    if sysinfo(2).contains("zh-cn")
        lang = "zh_CN";
    else
        lang = "en_US";
    end
end

end

