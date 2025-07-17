function is_safe = ValidatePluginSecurity(file)
%VALIDATEPLUGINSECURITY This function analysis code in file and report if
%there is dangerous behavior may destroy ReTiNA even MATLAB or computer.

arguments
    file    (1,1)   string
end

% TODO: 
% Defination of dangerous behaviour:  (should be external data base)
% (A) codes try to execute unsafe external library or destroy Java stack/heap
% (B) codes try to modify the properties constrained
% (C) codes try to control ReTiNA public UI components
% (D) codes try to run unacceptable compute resources with so much time/memory cost

is_safe = true;

end

