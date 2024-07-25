function info = GetRegionInfo()
%GETREGIONINFO This function get the region setting by using systeminfo
% at Windows platform or locale at unix platform
% Input:
% Output:
%   - info: struct, with platform and language setting

info = struct("plat",[],"lang",[],"mapping",struct());

if ispc()
    info.plat = "Windows";
    [~, sysinfo] = system("systeminfo");
    sysinfo = lower(string(sysinfo));
    if sysinfo.contains("zh-cn")
        info.lang = "zh_CN";
    elseif sysinfo.contains("en-us")
        info.lang = "en_US";
    elseif sysinfo.contains("ru-ru")
        info.lang = "ru_RU";
    elseif sysinfo.contains("fr-fr")
        info.lang = "fr_FR";
    elseif sysinfo.contains("es-es")
        info.lang = "es_ES";
    else
        info.lang = "en_US";
    end
elseif isunix()
    info.plat = "Unix";
    [~, sysinfo] = system("locale");
    sysinfo = lower(string(sysinfo).split("="));
    if sysinfo(2).contains("zh-cn")
        info.lang = "zh_CN";
    elseif sysinfo(2).contains("en-us")
        info.lang = "en_US";
    elseif sysinfo(2).contains("ru-ru")
        info.lang = "ru_RU";
    elseif sysinfo(2).contains("fr-fr")
        info.lang = "fr_FR";
    elseif sysinfo(2).contains("es-es")
        info.lang = "es_ES";
    else
        info.lang = "en_US";
    end
else
    info.plat = "Unknown";
    info.lang = "en_US";
    warning("Unsupported Operation System.");
end

end

