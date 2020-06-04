extract_CONN_atlas_FC(atlas,conn_dir,dataset)

% Extract functional connectivity matrices from an atlas within a CONN project
% Inputs:
% atlas (required)      : name of atlas within CONN project
% conn_dir (required)   : full path to Conn project directory
% dataset (required)    : name of dataset folder name

%% Settings
global globalDataDir; % this must be set within startup.m
datapath=[globalDataDir];

% Defaults
global globalAaronDataDir;
cd([globalAaronDataDir filesep 'ADHD_MIT']);

% make output directories
for i=1:length(sublist)
    cd([globalAaronDataDir filesep 'ADHD_MIT' filesep sublist{i}]);
    mkdir(['Extracted_ts_Conn']); 
end

for i=1:length(sublist)
    for j=1:1
    atlas_ts=[]; DMN_ts=[]; atlasNetworks_ts=[]; first_3min=[]; second_3min=[];
    clear data data_sessions xyz roi_files
    cd([globalAaronDataDir filesep 'ADHD_MIT' filesep 'conn_ADHD_MIT' filesep 'results/preprocessing']); 
    % find all subject files
    roi_files=dir(['*Condition000.mat']);
    load([roi_files(i).name]);
        for k=1:length(data) % find atlas ROIs
            if contains(names{k},atlas)
               atlas_ts=[atlas_ts data{k}]; 
            end
        end
        
    % extract FD (Jenkinson)
    cd([globalAaronDataDir filesep dataset filesep sublist{i}]);
    fd=dir(['*FDjenkinson*']);
    fd=load(fd.name); fd=fd.R;
    mean_FD=mean(fd(2:end)); mean_FD=mean_FD';
    
        % split into first vs. second 3-mins
    first_3min=atlas_ts(round(1:180/TR),:);
    second_3min=atlas_ts(round(180/TR+1:360)/TR,:);
    % save    
    cd([globalAaronDataDir filesep 'ADHD_MIT' filesep sublist{i} filesep 'Extracted_ts_Conn']);
    save('mean_FD','mean_FD');
    save([atlas '.txt'],'atlas_ts','-ascii');   
    save([atlas '_first_3min.txt'],'first_3min','-ascii'); 
    save([atlas '_second_3min.txt'],'second_3min','-ascii'); 
    display(['Done subject ' num2str(i) ' run ' num2str(j)]);
    end
end
