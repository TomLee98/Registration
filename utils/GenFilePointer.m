function [ma_ptr, ms_ptr] = GenFilePointer(mov, chl)
arguments
    mov (:,:,:,:,:) uint16
    chl (1,2) double {mustBeInRange(chl,1,3)} = [1, 2];
end

assert(chl(1)~=chl(2), "GenFilePointer:invalidChannelOrder", ...
    "The channel order is invalid.");

ma_ptr = [];
ms_ptr = [];

% tmp file folder
folder_unix = "/data/tmpdata/";
folder_pc = "tmpdata\";
folder_local = [];
filename.aligned = "mov_aligned_source.mat";
filename.signal = "mov_signal_source.mat";

% output file and handle the pointer
if isunix()
    [~, user_name] = system('id -u --name');
    user_name = string(user_name);
    aligned_file = folder_unix + user_name + "/" + filename.aligned;
    signal_file = folder_unix + user_name + "/" + filename.signal;
    folder_local = folder_unix + user_name + "/";
elseif ispc()
    user_name = string(getenv('username'));
    aligned_file = folder_pc + user_name + "\" + filename.aligned;
    signal_file = folder_pc + user_name + "\" + filename.signal;
    folder_local = folder_pc + user_name + "\";
end

if exist(aligned_file,"file") && exist(signal_file,"file")
    % ask for load the same name file?
    overlap_flag = input("There are memory mapping files, loading?(Y/N)","s");
    switch upper(overlap_flag)
        case "Y"
            % just loading and skip regenerating
            % establish file pointer
            disp("Loading mapping file...");
            ma_ptr = matfile(aligned_file,"Writable",true);
            ms_ptr = matfile(signal_file,"Writable",true);

            % flag for whether the data pointer loading into memory
            loading_flag = true;
        case "N"
            % goto the usual pipeline: saving file and reloading pointer
            % remove all exist temporary files
            delete(folder_local+"*.mat");

            loading_flag = false;
        otherwise
            error("unsupported operation.");
    end
else
    loading_flag = false;
end

if loading_flag == false
    disp("Memory mapping...");
    % save data to disk for big data support
    % ============= PROCESSING THE ALIGNED DATA ==============
    mov_aligned = squeeze(uint16(mov(:,:,chl(1),:,:)));
    if ~exist(folder_local, "dir")
        mkdir(folder_local);
    end
    % nocompression for faster saving/loading but more disk
    % space allocated
    save(aligned_file,"mov_aligned","-v7.3","-nocompression");
    clearvars mov_aligned;

    % ============= PROCESSING THE SIGNAL DATA ===============
    mov_signal = squeeze(uint16(mov(:,:,chl(2),:,:)));
    save(signal_file,"mov_signal","-v7.3","-nocompression");
    clearvars mov_signal;

    % create matfile object for dynamic matrix processing
    ma_ptr = matfile(aligned_file,"Writable",true);
    ms_ptr = matfile(signal_file,"Writable",true);
end

disp("Memory reloading succeed.");
end
