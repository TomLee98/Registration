function mt = modifySSD(ssd_data, funstr)
arguments
    ssd_data table
    funstr (1,1) string {mustBeMember(funstr, ["FlowChecker", "Others"])}...
        = "Others";
end
% This function only takes the rows without water
chemical_names = string(ssd_data.Var1);
if funstr == "Others"
    stay_rows = chemical_names.contains("uM")&(~chemical_names.contains("water"));
else
    % flow checker need consider the water flow speed 
    stay_rows = chemical_names.contains("uM");
end
% hold the 1st and 2nd rows on
stay_rows(1:2) = true;

% check the data format
valve_c = string(ssd_data.Var3(1));
if valve_c == "V" || valve_c == ""
    % This is version 2 file format
    % col2: concentration
    % col4: pressure port
    mt = ssd_data(stay_rows, [1,3,5:size(ssd_data, 2)]);
else
    mt = ssd_data(stay_rows,:);
end
end

