function info = getCpuInfo()
% code source:
% https://blog.csdn.net/qingfengxd1/article/details/124204794

if isunix % 判断matlab的版本 UNIX
    if ismac % 判断matlab的版本 Mac
        info = cpuInfoMac();
    else
        info = cpuInfoUnix();
    end
else
    info = cpuInfoWindows();
end
end


%-------------------------------------------------------------------------%
function info = cpuInfoWindows()
sysInfo = callWMIC('cpu');
osInfo = callWMIC('os');

info = struct( ...
    'Name', sysInfo.Name, ...
    'Clock', [sysInfo.MaxClockSpeed,' MHz'], ...
    'Cache', [sysInfo.L2CacheSize,' KB'], ...
    'NumProcessors', str2double( sysInfo.NumberOfCores ), ...
    'OSType', 'Windows', ...
    'OSVersion', osInfo.Caption ); %输出的字符串结构
end

%---------------------111111----------------------------------------------------%
function info = callWMIC( alias )
% 呼叫ms - dos WMIC(Windows管理)命令
olddir = pwd();  %显示当前文件夹
cd(tempdir);  % 计算机目录
sysinfo = evalc(sprintf('!wmic %s get /value',alias)); %得到计算机的参数
cd(olddir);
fields = textscan(sysinfo, '%s', 'Delimiter', '\n');
fields = fields{1};%生成cell 数据处理
fields(cellfun('isempty', fields)) = [];
% 每一行有“字段=值”,所以分开
values = cell( size( fields ) );
for ff=1:numel( fields )
    idx = find( fields{ff}=='=', 1, 'first' ); %找到一个first字符串
    if ~isempty( idx ) && idx>1
        values{ff} = strtrim( fields{ff}(idx+1:end) ); %字符串分割
        fields{ff} = strtrim( fields{ff}(1:idx-1) );
    end
end

% 删除任何重复（仅适用于双插槽电脑，我们将假设所有插座有相同的处理器）。
numResults = sum( strcmpi( fields, fields{1} ) );
if numResults>1
    % 计算核数
    numCoresEntries = find( strcmpi( fields, 'NumberOfCores' ) );
    if ~isempty( numCoresEntries )
        cores = cellfun( @str2double, values(numCoresEntries) );
        values(numCoresEntries) = {num2str( sum( cores ) )};
    end
    % 去掉重复结果
    [fields,idx] = unique(fields,'first');
    values = values(idx);
end

% 数据转换
info = cell2struct( values, fields );
end

%-------------------22222222------------------------------------------------------%
function info = cpuInfoMac()  %Mac和window类似
machdep = callSysCtl('machdep.cpu');
hw = callSysCtl('hw');
info = struct( ...
    'Name', machdep.brand_string, ...
    'Clock', [num2str(str2double(hw.cpufrequency_max)/1e6),' MHz'], ...
    'Cache', [machdep.cache.size,' KB'], ...
    'NumProcessors', str2double( machdep.core_count ), ...
    'OSType', 'Mac OS/X', ...
    'OSVersion', getOSXVersion() );
end
%-------------------------------------------------------------------------%
function info = callSysCtl( namespace )
infostr = evalc(sprintf('!sysctl -a %s', namespace));
% Remove the prefix
infostr = strrep(infostr, [namespace,'.'], ''); %字符串
% Now break into a structure
infostr = textscan(infostr, '%s', 'delimiter', '\n');
infostr = infostr{1};
info = struct();
for ii=1:numel(infostr)
    colonIdx = find( infostr{ii}==':', 1, 'first');
    if isempty( colonIdx ) || colonIdx==1 || colonIdx==length(infostr{ii})
        continue
    end
    prefix = infostr{ii}(1:colonIdx-1);
    value = strtrim(infostr{ii}(colonIdx+1:end));
    while ismember('.', prefix)
        dotIndex = find(prefix=='.', 1, 'last');
        suffix = prefix(dotIndex+1:end);
        prefix = prefix(1:dotIndex-1);
        value = struct(suffix, value);
    end
    info.(prefix) = value;

end
end
%-------------------3333333333------------------------------------------------------%
function vernum = getOSXVersion()%版本号
% 提取系统软件版本的操作系统版本号输出
ver = evalc('system(''sw_vers'')');%版本 执行matlab字符串 软件版本信息
vernum = regexp(ver, 'ProductVersion:\s([1234567890.]*)', 'tokens', 'once');%对字符串进行查找替换
vernum = strtrim(vernum{1});
end
%-------------------------------------------------------------------------%
function info = cpuInfoUnix()
txt = readCPUInfo();
cpuinfo = parseCPUInfoText(txt);

osinfo = parseOSInfoText();

% Merge the structures
info = cell2struct([struct2cell(cpuinfo);struct2cell(osinfo)], ...
    [fieldnames(cpuinfo);fieldnames(osinfo)]);
end
%-------------------------------------------------------------------------%
function info = parseCPUInfoText(txt)
% Now parse the fields
lookup = {
    'model name', 'Name'
    'cpu Mhz', 'Clock'
    'cpu cores', 'NumProcessors'
    'cache size', 'Cache'
    'physical id', 'NumCPUs'
    };
info = struct( ...
    'Name', {''}, ...
    'Clock', {''}, ...
    'Cache', {''},...
    'NumCPUs',{''});
for ii=1:numel(txt)
    if isempty(txt{ii})
        continue;
    end
    % Look for the colon that separates the property name from the value
    colon = find(txt{ii}==':', 1, 'first');
    if isempty(colon) || colon==1 || colon==length(txt{ii})
        continue;
    end
    fieldName = strtrim(txt{ii}(1:colon-1));
    fieldValue = strtrim(txt{ii}(colon+1:end));
    if isempty(fieldName) || isempty(fieldValue)
        continue;
    end

    % Is it one of the fields we're interested in?
    idx = find(strcmpi(lookup(:,1), fieldName));
    if ~isempty(idx)
        newName = lookup{idx,2};
        info.(newName) = fieldValue;
    end
end
info.NumCPUs = str2double(info.NumCPUs) + 1;

% Convert clock speed
info.Clock = [info.Clock, ' MHz'];

% Convert num cores
info.NumProcessors = str2double(info.NumProcessors)*info.NumCPUs;
end
%-----------------444444-------------------------------------------------%
function info = parseOSInfoText()
info = struct( ...
    'OSType', 'Linux');
end
%-------------------------------------------------------------------------%
function txt = readCPUInfo()

fid = fopen( '/proc/cpuinfo', 'rt' );
if fid<0
    error( 'cpuinfo:BadPROCCPUInfo', 'Could not open /proc/cpuinfo for reading' );
end
onCleanup(@()fclose(fid));

txt = textscan( fid, '%s', 'Delimiter', '\n' );
txt = txt{1};
end
%-------------------------------------------------------------------------%