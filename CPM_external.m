%function [r_p_all,data]=external_MW_CPM_rest(all_mats,all_behav,subs,...
%    partial_var,corr_method,FD,plotting,group_id,univariate,MW_saCPM,train_suffix)

function [r_p_all]=CPM_external(all_mats,all_behav,mdl,pos_mask,neg_mask,part_var,motion_var,r_method)

% Test CPM in external dataset 
% INPUTS:
% all_mats (required)   : ROI x ROI x trials FC matrix (or single vector for one ROI/edge)
% all_behav (required)  : behavioral score vector
% mdl (required)        : Coefficient fits for linear model (two values)
% pos_mask (required)   : mask for model's significant positive features
% neg_mask (required)   : mask for model's significant negative features
% part_var (optional)   : partial corr variable (leave blank if not using)
% motion_var (optional) : head motion as FD (if included, removes subjects with FD>0.15)
% r_method (optional)   : correlation method (1 = Pearson (default); 2 = spearman;
%                       3 = robust regress; 4 = Pearson partial using part_var;
%                       5 = Spearman partial using part_var    

% OUTPUTS:
% 

%% Settings
FD_thr=.15; % cutoff for remoing subjects based on FD

%% Defaults
global globalDataDir;
load('cdcol.mat');
datapath=[globalDataDir filesep 'MW_gradCPT']; % location of training data
%cd([datapath]);

if nargin<8 || isempty(group_id)
   group_id=[]; 
end
if nargin<5 || isempty(corr_method)
   corr_method=1; 
end

if nargin<4 || isempty(partial_var)
 partial_var=[];
end
if nargin<6 || isempty(FD)
 FD=[];
end
if nargin<7 || isempty(plotting)
 plotting=1;
end
if nargin<10 || isempty(MW_saCPM)
MW_saCPM=1;
end
if nargin<9 || isempty(univariate)
univariate=[];
end
if nargin<11 || isempty(train_suffix)
train_suffix='HC';
end

% get fit parameters, positive and negative masks from training data
%cd([datapath '/Group/Results']);
if MW_saCPM==1
load([datapath '/Group/Results/fit_posneg_' train_suffix '.mat']);
load([datapath '/Group/Results/fit_pos_' train_suffix '.mat']);
load([datapath '/Group/Results/fit_neg_' train_suffix '.mat']);
load([datapath '/Group/Results/pos_mask_' train_suffix '.mat']);
load([datapath '/Group/Results/neg_mask_' train_suffix '.mat']);
elseif MW_saCPM==2
load([datapath '/Group/Results/fit_posneg_ADHD.mat']);
load([datapath '/Group/Results/fit_pos_ADHD.mat']);
load([datapath '/Group/Results/fit_neg_ADHD.mat']);
load([datapath '/Group/Results/pos_mask_ADHD.mat']);
load([datapath '/Group/Results/neg_mask_ADHD.mat']);   
elseif MW_saCPM==3
load([datapath '/Group/Results/high_attention_mask_saCPM.mat']);
load([datapath '/Group/Results/low_attention_mask_saCPM.mat']);   
load([datapath '/Group/Results/fit_posneg_saCPM.mat']);     
end

% set parameters
no_sub=length(subs);

if ndims(all_mats)==3
    no_node=size(all_mats,1);
    sub_trials=1:size(all_mats,3);
end
no_trials_all=length(all_behav);
no_trials=length(all_behav)/no_sub;

behav_pred_pos=zeros(no_trials,1);
behav_pred_neg=zeros(no_trials,1);

% remove high-motion subjects
if ~isempty(FD) 
    rm_subs=find(FD>FD_thr);
    if corr_method~=3
    display(['removing ' num2str(length(rm_subs)) ' subjects due to high motion']);
    all_behav(rm_subs)=NaN;
    partial_var(rm_subs)=NaN;
    end
end

curr_sub_trials=1;
for leftout=1:no_sub
    % extract one subject from matrices and behavior   
    % run model on test subject
    if ndims(all_mats)==3
    test_mat=all_mats(:,:,leftout);
    test_sumpos(leftout)=nansum(nansum(test_mat.*pos_mask))/2;
    test_sumneg(leftout)=nansum(nansum(test_mat.*neg_mask))/2;
    test_sum_posneg(leftout)=squeeze(test_sumpos(leftout))-squeeze(test_sumneg(leftout));
    
    behav_pred_posneg(leftout)=fit_posneg(1)*test_sum_posneg(leftout) + fit_posneg(2);
    if MW_saCPM<3
    behav_pred_pos(leftout)=fit_pos(1)*test_sumpos(leftout) + fit_pos(2);
    behav_pred_neg(leftout)=fit_neg(1)*test_sumneg(leftout) + fit_neg(2);
    end
    end
    curr_sub_trials=curr_sub_trials+no_trials;
end
    
% compare predicted and observed scores
if corr_method<3 || corr_method==4
    if MW_saCPM<3
[R_pos,P_pos]=corr(behav_pred_pos',all_behav,'rows','pairwise');
[R_neg,P_neg]=corr(behav_pred_neg',all_behav,'rows','pairwise');
    end
[R_posneg,P_posneg]=corr(behav_pred_posneg',all_behav,'rows','pairwise');
[spearman_R_posneg,spearman_P_posneg]=corr(behav_pred_posneg',all_behav,'type','Spearman','rows','pairwise');
if MW_saCPM<3
[spearman_R_pos,spearman_P_pos]=corr(behav_pred_pos',all_behav,'type','Spearman','rows','pairwise');
[spearman_R_neg,spearman_P_neg]=corr(behav_pred_neg',all_behav,'type','Spearman','rows','pairwise');
end

if ~isempty(partial_var)
    if MW_saCPM<3
   [partial_R_pos,partial_P_pos]=partialcorr(behav_pred_pos',all_behav,partial_var,'rows','pairwise');
   [partial_R_neg,partial_P_neg]=partialcorr(behav_pred_neg',all_behav,partial_var,'rows','pairwise');
    end
   [partial_R_posneg,partial_P_posneg]=partialcorr(behav_pred_posneg',all_behav,partial_var,'rows','pairwise'); 
   if MW_saCPM<3
   [partial_rho_pos,partial_rho_P_pos]=partialcorr(behav_pred_pos',all_behav,partial_var,'type','Spearman','rows','pairwise');
   [partial_rho_neg,partial_rho_P_neg]=partialcorr(behav_pred_neg',all_behav,partial_var,'type','Spearman','rows','pairwise');
   end
   [partial_rho_posneg,partial_rho_P_posneg]=partialcorr(behav_pred_posneg',all_behav,partial_var,'type','Spearman','rows','pairwise'); 
end
end

if corr_method==3 % if LME 
    % LME: posneg CPM vs behavior
    sub_code=cell2mat(subs);
    lmm_behav=all_behav;
    lmm_variables=[sub_code all_behav test_sum_posneg'];
    lmm_table=table(lmm_variables(:,1),lmm_variables(:,2),lmm_variables(:,3),...
    'VariableNames',{'subject','behav','FC'});
    lme=fitlme(lmm_table,'behav~FC+(1|subject)');
    lme_anova=anova(lme,'DFmethod','satterthwaite');
    
    % compare runs
    test_sum_posneg_runs=[];
    if ~isempty(FD)
        test_sum_posneg(rm_subs)=NaN;
    end
    for i=1:max(sub_code)
        test_sum_posneg_runs(i,:)=test_sum_posneg(find(sub_code==i));
    end
    test_sum_posneg_runs_norm=[];
    for i=1:size(test_sum_posneg_runs,1)
        test_sum_posneg_runs_norm(i,:)=zscore(test_sum_posneg_runs(i,:)); 
    end
    test_sum_posneg_norm=[];
    for i=1:size(test_sum_posneg_runs_norm,2)
        test_sum_posneg_norm=[test_sum_posneg_norm; test_sum_posneg_runs_norm(:,i)];
    end
    all_behav_runs=[];
    for i=1:max(sub_code)
        all_behav_runs(i,:)=all_behav(find(sub_code==i));
    end
    % LME: interaction between CPM posneg and run number
    sub_code=cell2mat(subs);
    lmm_cpm=test_sum_posneg';
    lmm_run=[];
    for i=1:size(test_sum_posneg_runs,2)
       lmm_run=[lmm_run; [i*ones(size(test_sum_posneg_runs,1),1)]];
    end
    cpm_lmm_variables=[sub_code lmm_run lmm_cpm];
    cpm_lmm_table=table(cpm_lmm_variables(:,1),cpm_lmm_variables(:,2),cpm_lmm_variables(:,3),...
    'VariableNames',{'subject','run','CPM'});
    cpm_lme=fitlme(cpm_lmm_table,'CPM~run+(1|subject)');
    cpm_lme_anova=anova(cpm_lme,'DFmethod','satterthwaite');
    
    % LME: interaction between behav and run number 
    lmm_behav=all_behav;
    behav_lmm_variables=[sub_code lmm_run lmm_behav];
    behav_lmm_table=table(behav_lmm_variables(:,1),behav_lmm_variables(:,2),behav_lmm_variables(:,3),...
    'VariableNames',{'subject','run','behav'});
    behav_lme=fitlme(behav_lmm_table,'behav~run+(1|subject)');
    behav_lme_anova=anova(behav_lme,'DFmethod','satterthwaite');
end

if isempty(univariate)==0
    [R_uni,P_uni]=corr(univariate,all_behav,'rows','pairwise');
    [spearman_R_uni,spearman_P_uni]=corr(univariate,all_behav,'type','Spearman','rows','pairwise');
    [partial_R_uni,partial_P_uni]=partialcorr(univariate,all_behav,partial_var,'rows','pairwise');
    [partial_rho_uni,partial_rho_P_uni]=partialcorr(univariate,all_behav,partial_var,'type','Spearman','rows','pairwise');
end

%% organize output (R and P values)
if corr_method==4
r_p_all=[R_posneg P_posneg spearman_R_posneg spearman_P_posneg partial_R_posneg partial_P_posneg partial_rho_posneg partial_rho_P_posneg];
data=[all_behav behav_pred_posneg' test_sum_posneg'];
elseif corr_method==1 || corr_method==2
r_p_all=[R_posneg P_posneg];
data=[all_behav behav_pred_posneg' test_sum_posneg'];
else
r_p_all=[];
data.all_behav_runs=all_behav_runs;
data.behav_lme_anova=behav_lme_anova;
data.test_sum_posneg_runs=test_sum_posneg_runs;
data.test_sum_posneg_runs_norm=test_sum_posneg_runs_norm;
end

%% ploting (call functions)
if isempty(univariate)==1
    if plotting==1
if corr_method==1 || corr_method==2 || corr_method==4
savepath=[datapath filesep 'Group' filesep 'Results' filesep 'figs' filesep];
if isempty(group_id)==0
    MW_CPM_plot_external_rest_subgroups(savepath,group_id,all_behav,partial_var,behav_pred_posneg,spearman_R_posneg,spearman_P_posneg,...
    partial_R_posneg,partial_P_posneg,partial_rho_posneg,partial_rho_P_posneg)

else
        MW_CPM_plot_external_rest(MW_saCPM,savepath,all_behav,partial_var,behav_pred_posneg,behav_pred_pos,...
            behav_pred_neg,spearman_R_posneg,spearman_P_posneg,...
    partial_R_posneg,partial_P_posneg,partial_rho_posneg,partial_rho_P_posneg,R_posneg,P_posneg,R_pos,P_pos,spearman_R_pos,...
    spearman_P_pos,partial_R_pos,partial_P_pos,partial_rho_pos,partial_rho_P_pos,R_neg,P_neg,spearman_R_neg,...
    spearman_P_neg,partial_rho_neg,partial_rho_P_neg)
end
end
if corr_method==3
 MW_CPM_plot_lmm_runs(behav_lme_anova,test_sum_posneg_runs,test_sum_posneg_runs_norm,all_behav_runs);   
end
elseif isempty(univariate==0)
    
end
end
