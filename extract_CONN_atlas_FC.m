function extract_CONN_atlas_FC(atlas,sublist,conn_dir,dataset)

% Extract functional connectivity matrices from an atlas within a CONN project
% INPUTS:
% atlas (required)      : name of atlas within CONN project (e.g.'Schaefer300')
% sublist (required)    : cell array with folders names of all subjects included in CONN project
% conn_dir (required)   : full path to Conn project directory
% dataset (required)    : name of dataset folder name
% OUTPUTS:
% within Extracted_ts_Conn subfolder for each subject: FC matrices
% within dataset/Group folder: saves .mat file with 3D matrix (node x node x subject) r; also saves vector of FD (head motion)

%% Settings
global globalDataDir; % this must be set within startup.m
datapath=[globalDataDir];

% make output directories
for i=1:length(sublist)
    mkdir([datapath filesep dataset filesep sublist{i} filesep 'Extracted_ts_Conn']); 
end
mkdir([datapath filesep dataset filesep 'Group']); 

% load conn project file
a=strfind(conn_dir,filesep);
pre_conn_dir=conn_dir(1:a(end));
post_conn_dir=conn_dir(a(end)+1:length(conn_dir));
load([pre_conn_dir post_conn_dir '.mat']);

% extract time series
for i=1:length(sublist)
    for j=1:1 % edit for >1 run?
    atlas_ts=[];
    clear data data_sessions xyz roi_files
    % find all subject files
    roi_files=dir([conn_dir filesep 'results' filesep 'preprocessing' filesep '*Condition000.mat']);
    load([conn_dir filesep 'results' filesep 'preprocessing' filesep roi_files(i).name]);
        for k=1:length(data) % find atlas ROIs
            if contains(names{k},atlas)
               atlas_ts=[atlas_ts data{k}]; 
            end
        end
        
    % extract FD (Jenkinson)
    b=[];
    FD_dir=CONN_x.Setup.functional{1,i}{1,1}{1,3}.fname; % edit for >1 run?   
    b=strfind(FD_dir,filesep);
    fd_dir=FD_dir(1:b(end));
    fd=dir([fd_dir  '*FDjenkinson*']);
    fd=load([fd_dir fd.name]); 
    fd=fd.R;
    mean_FD=mean(fd(2:end)); mean_FD=mean_FD';

    % save subject data    
    save([datapath filesep dataset filesep sublist{i} filesep 'Extracted_ts_Conn' filesep 'mean_FD'],'mean_FD');
    save([datapath filesep dataset filesep sublist{i} filesep 'Extracted_ts_Conn' filesep atlas '.txt'],'atlas_ts','-ascii');   
    display(['Done subject ' num2str(i) ' run ' num2str(j)]);
    end
end

% merge across subjects
display(['Merging subjects : see output at ' datapath '/' dataset '/Group']);
mean_FD=[];
for i=1:length(sublist)
    atlas_ts=[]; 
    atlas_ts=importdata([globalDataDir filesep dataset filesep sublist{i} '/Extracted_ts_Conn' filesep atlas '.txt']);
    atlas_all_mats(:,:,i)=fisherz(corrcoef(atlas_ts));
    load([globalDataDir filesep dataset filesep sublist{i} '/Extracted_ts_Conn' filesep 'mean_FD.mat']);
    mean_FD_all(i)=mean_FD;
end
mean_FD=mean_FD_all';

% save concatenated FC matrices and FD
save([datapath filesep dataset filesep 'Group' filesep atlas '_all_mats'],'atlas_all_mats');
save([datapath filesep dataset filesep 'Group' filesep 'mean_FD','mean_FD']);
