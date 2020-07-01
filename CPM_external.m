function [r_p_all,data]=CPM_external(all_mats,all_behav,cpm,part_var,motion_var)

% Test CPM in external dataset 
% written by Aaron Kucyi, Northeastern University
% INPUTS:
% all_mats (required)   : ROI x ROI x trials FC matrix (or single vector
%                       for one ROI/edge) from test dataset
% all_behav (required)  : behavioral score vector from test dataset
% cpm (required)        : structure containing coefficient fits (fit_posneg), pos_mask and neg_mask 
%                       from training dataset
% part_var (optional)   : partial corr variable (leave blank if not using)
% motion_var (optional) : head motion as FD (if included, removes subjects with FD>0.15)

% OUTPUTS:
% r_p_all = r (Pearson), p (Pearson), rho (Spearman), p (Spearman)
%           [plus repeat for partial correlations if selected]
% data: predicted behav (column 1), observed behav (column 2), network strength (column 3)

%% Settings
FD_thr=.15; % cutoff for remoing subjects based on FD

%% Defaults
global globalDataDir;
if nargin<4 || isempty(part_var)
   part_var=[]; 
end
if nargin<5 || isempty(motion_var)
   motion_var=[]; 
end

% set parameters
no_sub=length(all_behav);
if ndims(all_mats)==3
    no_node=size(all_mats,1);
    sub_trials=1:size(all_mats,3);
end

behav_pred_pos=zeros(no_sub,1);
behav_pred_neg=zeros(no_sub,1);

% remove high-motion subjects
if ~isempty(motion_var) 
    rm_subs=find(motion_var>FD_thr);
    display(['removing ' num2str(length(rm_subs)) ' subjects due to high motion']);
    all_behav(rm_subs)=NaN;
    if ~isempty(part_var)
    part_var(rm_subs)=NaN;
    end
end

% loop through test subjects
curr_sub_trials=1;
for leftout=1:no_sub
    if ndims(all_mats)==3
    test_mat=all_mats(:,:,leftout);
    test_sumpos(leftout)=nansum(nansum(test_mat.*cpm.pos_mask))/2;
    test_sumneg(leftout)=nansum(nansum(test_mat.*cpm.neg_mask))/2;
    test_sum_posneg(leftout)=squeeze(test_sumpos(leftout))-squeeze(test_sumneg(leftout));
    behav_pred_posneg(leftout)=cpm.fit_posneg(1)*test_sum_posneg(leftout) + cpm.fit_posneg(2);
    behav_pred_pos(leftout)=cpm.fit_pos(1)*test_sumpos(leftout) + cpm.fit_pos(2);
    behav_pred_neg(leftout)=cpm.fit_neg(1)*test_sumneg(leftout) + cpm.fit_neg(2); 
    end
    curr_sub_trials=curr_sub_trials+no_sub;
end
    
% compare predicted and observed scores
[R_pos,P_pos]=corr(behav_pred_pos,all_behav,'rows','pairwise');
[R_neg,P_neg]=corr(behav_pred_neg,all_behav,'rows','pairwise');
[R_posneg,P_posneg]=corr(behav_pred_posneg',all_behav,'rows','pairwise');
[spearman_R_posneg,spearman_P_posneg]=corr(behav_pred_posneg',all_behav,'type','Spearman','rows','pairwise');
[spearman_R_pos,spearman_P_pos]=corr(behav_pred_pos,all_behav,'type','Spearman','rows','pairwise');
[spearman_R_neg,spearman_P_neg]=corr(behav_pred_neg,all_behav,'type','Spearman','rows','pairwise');

if ~isempty(part_var)
   [partial_R_pos,partial_P_pos]=partialcorr(behav_pred_pos,all_behav,part_var,'rows','pairwise');
   [partial_R_neg,partial_P_neg]=partialcorr(behav_pred_neg,all_behav,part_var,'rows','pairwise');
   [partial_R_posneg,partial_P_posneg]=partialcorr(behav_pred_posneg',all_behav,part_var,'rows','pairwise'); 
   [partial_rho_pos,partial_rho_P_pos]=partialcorr(behav_pred_pos,all_behav,part_var,'type','Spearman','rows','pairwise');
   [partial_rho_neg,partial_rho_P_neg]=partialcorr(behav_pred_neg,all_behav,part_var,'type','Spearman','rows','pairwise');
   [partial_rho_posneg,partial_rho_P_posneg]=partialcorr(behav_pred_posneg',all_behav,part_var,'type','Spearman','rows','pairwise'); 
end

%% organize output (R and P values, predicted vs observed vs network strength)
if ~isempty(part_var)
r_p_all=[R_posneg P_posneg spearman_R_posneg spearman_P_posneg partial_R_posneg partial_P_posneg partial_rho_posneg partial_rho_P_posneg];
data=[behav_pred_posneg' all_behav test_sum_posneg'];
else
r_p_all=[R_posneg P_posneg spearman_R_posneg spearman_P_posneg];
data=[behav_pred_posneg' all_behav test_sum_posneg'];
end

