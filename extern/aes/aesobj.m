function obj = aesobj()
%AESOBJ This function use a fix code as key and generate the aes object
% The file will be generate as pcode for hidding the details

% key components:
% mother board_cpu_hard drive_network interface controller

key = genkey("MB|CPU|HDD/SSD|NIC");

obj = AES(string(key), "SHA-1");
end

function key = genkey(comp, delimiter)
arguments
    comp (1,1) string = "MB|CPU|HDD/SSD";
    delimiter (1,1) string = "_";
end

comp = comp.split("|");
key = strings(numel(comp), 1);

for k = 1:numel(comp)
    switch comp(k)
        case "MB"
            cmd = 'wmic baseboard get serialnumber';
        case "CPU"
            cmd = 'wmic cpu get processorid';
        case "HDD/SSD"
            cmd = 'wmic diskdrive get serialnumber';
        case "NIC"
            cmd = 'wmic nicconfig get macaddress';
    end
    [~, result] = system(cmd);
    fields = textscan( result, '%s', 'Delimiter', '\n' );
    fields = strtrim(fields{1});
    fields = string(fields(2:end));
    fields(fields == "") = [];

    key(k) = fields.join("_");
end

key = key.join(delimiter);
end