function status = CompileLanguages(lang)
%GENERATEMULTILANGSUPPORT This function generate multiple languages
%supports, include "zh-CN","en-US","ru-RU","fr-FR","es-ES"
% Input:
%   - lang: string, can be "zh-CN","en-US","ru-RU","fr-FR","es-ES"
% Output:
%   - status: generation status

arguments
    lang (1,1) string ...
        {mustBeMember(lang,["zh_CN","en_US","all"])} = "en_US";
end

% load the identifier-language pair file
if ispc()
    p = readtable("lang\words_dict.xlsx","TextType","string",...
        "FileType","spreadsheet","ExpectedNumVariables",3);
elseif isunix()
    p = readtable("lang/words_dict.xlsx","TextType","string",...
        "FileType","spreadsheet","ExpectedNumVariables",3);
end

chmap = struct();

if lang == "all"
    langs = ["zh_CN", "en_US"];
else
    langs = lang;
end

for lang = langs
    for k = 1:size(p,1)
        chmap.(p.ID(k)) = p.(lang)(k);
    end
    if ispc()
        writestruct(chmap,"lang\"+lang+".xml","FileType","xml");
    elseif isunix()
        writestruct(chmap,"lang/"+lang+".xml","FileType","xml");
    end
end

status = 0;
end