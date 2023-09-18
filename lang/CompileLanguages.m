function status = CompileLanguages(lang)
%GENERATEMULTILANGSUPPORT This function generate multiple languages
%supports, include "zh-CN","en-US","ru-RU","fr-FR","es-ES"
% Input:
%   - lang: string, can be "zh-CN","en-US","ru-RU","fr-FR","es-ES"
% Output:
%   - status: generation status

arguments
    lang (1,1) string ...
        {mustBeMember(lang,["zh_CN","en_US","ru_RU","fr_FR","es_ES"])} = "en_US";
end

% load the identifier-language pair file
p = readtable("lang\words_dict.xlsx","TextType","string",...
    "FileType","spreadsheet","ExpectedNumVariables",6);

mapping = struct();

for k = 1:size(p,1)
     mapping.(p.ID(k)) = p.(lang)(k);
end

writestruct(mapping,"lang\"+lang+".xml","FileType","xml");

status = 0;
end

