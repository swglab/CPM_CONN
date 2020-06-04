function [R_posneg,P_posneg,behav_pred_posneg]=CPM_internal(all_mats,all_behav,dataset,r_method,pthresh,part_var,motion_var,outname,train_mode)

% Connectome-based predictive modeling: internal (within-dataset) validation 
% written by Aaron Kucyi, Northeastern University
% adapted from Shen et al. (2017 Nature Protocols)
% INPUTS:
% all_mats (required)   : ROI x ROI x trials FC matrix (or single vector for one ROI/edge)
% all_behav (required)  : behavioral score vector
% dataset (required)    : name of dataset folder name
% r_method (optional)   : correlation method (1 = Pearson (default); 2 = spearman;
%                       3 = robust regress; 4 = Pearson partial using part_var;
%                       5 = Spearman partial using part_var          
% pthresh (optional)    : p threshold for feature selection (default = 0.01)
% part_var (optional)   : partial corr variable (leave blank if not using)
% motion_var (optional) : head motion as FD (if included, removes subjects with FD>0.15)
% outname (optional)    : name for output files (default = 'test')
% train_mode (optional) : 1 = do cross-validation (default); 2 = use all subjects for
%                       training in a single fold and save linear regression parameters
% OUTPUTS:
% R_posneg              : r for predicted vs observed behav (pos-neg model)
% P_posneg              : p for predicted vs observed behav (pos-neg model)
% behav_pred_posneg     : vector of predicted behavior
% positive and negative masks are saved in cpm_results within dataset folder

%% Settings
FD_thr=.15; % cutoff for removing subjects based on FD
global globalDataDir; % this must be set within startup.m
datapath=[globalDataDir];

%% Defaults
if nargin<4 || isempty(r_method)
    r_method=1;
end
if nargin<5 || isempty(pthresh)
    pthresh=0.01;
end
if nargin<6 || isempty(part_var)
    part_var=[];
end
if nargin<7 || isempty(motion_var)
    motion_var=[];
end
if nargin<8 || isempty(outname)
    outname='test';
end
if nargin<9 || isempty(train_mode)
    train_mode=1;
end

%% Defaults
% set parameters
if train_mode==1
no_sub=length(all_behav);
elseif train_mode==2
no_sub=1; % set to one loop if using all subjects for training
end

if ndims(all_mats)==3
    no_node=size(all_mats,1);
    sub_trials=1:size(all_mats,3);
else % for single ROI
    no_node=1;
    sub_trials=1:length(all_mats);
end
behav_pred_pos=zeros(no_sub,1);
behav_pred_neg=zeros(no_sub,1);

%% remove high-motion subjects
if ~isempty(motion_var) 
    rm_subs=find(motion_var>FD_thr);
    if r_method~=3
    display(['removing ' num2str(length(rm_subs)) ' subjects due to high motion']);
    all_behav(rm_subs)=NaN;
    partial_var(rm_subs)=NaN;
    end
end

%% Run cross-validation folds
for leftout=1:no_sub
    if train_mode==1 % leave one subject out
        display(['Running fold ' num2str(leftout)]);   
        % leave one subject from FC matrices and behavior   
        train_mats = all_mats;
        if ndims(all_mats)==3
        train_mats(:,:,leftout)=[];
        train_vcts=reshape(train_mats,[],size(train_mats,3));
        else % for single ROI/edge
        train_vcts=train_mats;
        train_vcts(leftout)=[];
        end
        train_behav=all_behav;
        train_behav(leftout)=[];
    elseif train_mode==2
        display(['Using all subjects to estimate model parameters']);
        train_mats = all_mats;
        train_vcts=reshape(train_mats,[],size(train_mats,3));
        train_behav=all_behav;
    end
    
    % correlate all edges with behavior (in training data)
    if r_method==1 || r_method==4 || r_method==5 
        if ndims(all_mats)==3
    [r_mat,p_mat]=corr(train_vcts',train_behav,'rows','pairwise');
    r_mat=reshape(r_mat,no_node,no_node);
    p_mat=reshape(p_mat,no_node,no_node);
        else % for single ROI/edge
            [r_mat,p_mat]=corr(train_vcts,train_behav,'rows','pairwise');
        end
    elseif r_method==2
        if ndims(all_mats)==3
    [r_mat,p_mat]=corr(train_vcts',train_behav,'type','Spearman','rows','pairwise');
    r_mat=reshape(r_mat,no_node,no_node);
    p_mat=reshape(p_mat,no_node,no_node);
        else % for single ROI/edge
        [r_mat,p_mat]=corr(train_vcts,train_behav,'type','Spearman','rows','pairwise');    
        end
    elseif r_method==3 % robust regression
        edge_no=size(train_vcts,1);
        r_mat=zeros(1,edge_no);
        p_mat=zeros(1,edge_no);
        for edge_i = 1:edge_no
           [~,stats]=robustfit(train_vcts(edge_i,:)',train_behav); 
           cur_t = stats.t(2);
           r_mat(edge_i) = sign(cur_t)*sqrt(cur_t^2/no_trials_all-1-2+cur_t^2);
           p_mat(edge_i) = 2*tcdf(cur_t, no_trials_all-1-2); % two-tailed
        end
        r_mat=reshape(r_mat,no_node,no_node);
        p_mat=reshape(p_mat,no_node,no_node);
    end   
    
    % set threshold and define masks (in training data)
    if ndims(all_mats==3)
    pos_mask=zeros(no_node,no_node);
    neg_mask=zeros(no_node,no_node);
    
    pos_edges=find(r_mat>0 & p_mat < pthresh);
    neg_edges=find(r_mat<0 & p_mat < pthresh);
    
    pos_mask(pos_edges)=1;
    neg_mask(neg_edges)=1;
    
    pos_mask_allfolds(:,:,leftout)=pos_mask;
    neg_mask_allfolds(:,:,leftout)=neg_mask;
    end
    
    % get sum of all edges in train subs
    if ndims(all_mats==3)
        if train_mode==1
        train_sumpos=zeros(size(all_mats,3)-1,1); % -1 is for leave one out; edit for k-fold
        train_sumneg=zeros(size(all_mats,3)-1,1);
        train_sumposneg=zeros(size(all_mats,3)-1,1);
        elseif train_mode==2 
        train_sumpos=zeros(size(all_mats,3),1);
        train_sumneg=zeros(size(all_mats,3),1);
        train_sumposneg=zeros(size(all_mats,3),1);     
        end
    
        for ss=1:size(train_sumpos)
        train_sumpos(ss) = nansum(nansum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = nansum(nansum(train_mats(:,:,ss).*neg_mask))/2;  
        end
        train_sumposneg = train_sumpos-train_sumneg;
    end
    
    % build model on train subjects
    if ndims(all_mats)==3
    fit_pos=polyfit(train_sumpos(~isnan(train_behav)),train_behav(~isnan(train_behav)),1);
    fit_neg=polyfit(train_sumneg(~isnan(train_behav)),train_behav(~isnan(train_behav)),1);
    fit_posneg=polyfit(train_sumposneg(~isnan(train_behav)),train_behav(~isnan(train_behav)),1);
    else % for single ROI/edge and positive corr only
    fit_pos=polyfit(train_vcts(~isnan(train_vcts)),train_behav(~isnan(train_behav)),1);
    end
    
    % run model prediction on test subject
    if train_mode==1
    if ndims(all_mats)==3
    test_mat=all_mats(:,:,leftout);
    test_sumpos(leftout)=nansum(nansum(test_mat.*pos_mask))/2;
    test_sumneg(leftout)=nansum(nansum(test_mat.*neg_mask))/2;
    test_sum_posneg(leftout)=squeeze(test_sumpos(leftout))-squeeze(test_sumneg(leftout)); 
    behav_pred_pos(leftout)=fit_pos(1)*test_sumpos(leftout) + fit_pos(2); 
    behav_pred_neg(leftout)=fit_neg(1)*test_sumneg(leftout) + fit_neg(2);
    behav_pred_posneg(leftout,:)=fit_posneg(1)*test_sum_posneg(leftout) + fit_posneg(2);
    else
    test_mat=all_mats(leftout);
    for j=1:length(leftout)
    behav_pred_pos(leftout,:)=fit_pos(1)*test_mat + fit_pos(2); 
    end
    end
    end
end

% compare predicted vs. observed scores across subjects
if train_mode==1
behav_pred_pos(isnan(all_behav))=NaN;
behav_pred_neg(isnan(all_behav))=NaN;
behav_pred_posneg(isnan(all_behav))=NaN;

    if r_method==1 || r_method==2 || r_method==3
      if ndims(all_mats)==3  
        [R_pos, P_pos] = corr(behav_pred_pos,all_behav,'rows','pairwise');
        [R_neg, P_neg] = corr(behav_pred_neg,all_behav,'rows','pairwise');
        [R_posneg, P_posneg] = corr(behav_pred_posneg,all_behav,'rows','pairwise');
      else % for single ROI and positive corr
    [R_posneg(leftout), P_posneg(leftout)] = corr(behav_pred_pos,behav_test,'rows','pairwise');    
        end
    elseif r_method==4 % partial Pearson
    [R_pos, P_pos] = partialcorr(behav_pred_pos,all_behav,part_var,'rows','pairwise');    
    [R_neg, P_neg] = partialcorr(behav_pred_neg,all_behav,part_var,'rows','pairwise');
    [R_posneg, P_posneg] = partialcorr(behav_pred_posneg,all_behav,part_var,'rows','pairwise');
    elseif r_method==5 % partial Spearman
    [R_pos, P_pos] = partialcorr(behav_pred_pos,all_behav,part_var,'type','Spearman','rows','pairwise');    
    [R_neg, P_neg] = partialcorr(behav_pred_neg,all_behav,part_var,'type','Spearman','rows','pairwise');
    [R_posneg, P_posneg] = partialcorr(behav_pred_posneg,all_behav,part_var,'type','Spearman','rows','pairwise');
    end
end

% save linear regression parameters (when using all subjects for training),
% positive mask and negative mask
mkdir([datapath filesep dataset filesep 'cpm_results']);
savepath=[datapath filesep dataset filesep 'cpm_results'];
if train_mode==2
save([savepath filesep 'fit_posneg_' outname],'fit_posneg');
save([savepath filesep 'fit_pos_' outname],'fit_pos');
save([savepath filesep 'fit_neg_' outname],'fit_neg');
save([savepath filesep 'pos_mask_' outname],'pos_mask'); 
save([savepath filesep 'neg_mask_' outname],'neg_mask');
save([savepath filesep 'pos_mask_' outname '.txt'],'pos_mask','-ascii');
save([savepath filesep 'neg_mask_' outname '.txt'],'neg_mask','-ascii');
R_posneg=[]; P_posneg=[]; behav_pred_pos=[];
end

% find positive and negative edges with overlap across 90% of all folds, save
if ndims(all_mats)==3 && train_mode==1
pos_mask_overall=NaN(size(pos_mask_allfolds,1),size(pos_mask_allfolds,2));
for i=1:size(pos_mask_allfolds,2)
    for j=1:size(pos_mask_allfolds,2)
        if sum(pos_mask_allfolds(i,j,:))/size(pos_mask_allfolds,3)>.9
        pos_mask_overall(i,j)=1;
        else
        pos_mask_overall(i,j)=0;
        end
    end
end
save([savepath filesep 'pos_mask_90percFolds_' outname '.txt'],'pos_mask_overall','-ascii');
save([savepath filesep 'pos_mask_90percFolds' outname],'pos_mask_allfolds');

neg_mask_overall=NaN(size(neg_mask_allfolds,1),size(neg_mask_allfolds,2));
for i=1:size(neg_mask_allfolds,2)
    for j=1:size(neg_mask_allfolds,2)
        if sum(neg_mask_allfolds(i,j,:))/size(neg_mask_allfolds,3)>.9
        neg_mask_overall(i,j)=1;
        else
        neg_mask_overall(i,j)=0;
        end
    end
end
save([savepath filesep 'neg_mask_90percFolds_' outname '.txt'],'neg_mask_overall','-ascii');
save([savepath filesep 'neg_mask_90percFolds_' outname],'neg_mask_allfolds');
end

