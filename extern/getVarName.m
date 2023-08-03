function [str_varName] = getVarName(var) %#ok<INUSD>
%GETVARNAME This function output var name as string
str_varName = sprintf('%s',inputname(1));
end

