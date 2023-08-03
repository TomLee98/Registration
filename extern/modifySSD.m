function mt = modifySSD(ssd_data)
chemical_names = string(ssd_data.Var1);
stay_rows = chemical_names.contains("uM");
% hold the 1st and 2nd rows on
stay_rows(1:2) = true;

% check the data format
valve_c = string(ssd_data.Var3(1));
if valve_c == "V"
    % col2: concentration
    % col4: pressure port
    mt = ssd_data(stay_rows, [1,3,5:size(ssd_data, 2)]);
else
    mt = ssd_data(stay_rows,:);
end
end

