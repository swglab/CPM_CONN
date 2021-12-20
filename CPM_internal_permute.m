function [p_all,R,R_permute_all]=CPM_internal_permute(all_mats,all_behav,dataset,...
    kfolds,r_method,pthresh,part_var,motion_var,outname,no_iter)

% Permutation test for CPM_internal.m
% written by Aaron Kucyi, Northeastern University
% adapted from Shen et al. (2017 Nature Protocols)
% INPUTS:
% all_mats (required)   : ROI x ROI x trials FC matrix (or single vector for one ROI/edge)
% all_behav (required)  : behavioral score vector
% dataset (required)    : name of dataset folder name
% kfolds (optional)     : number of cross-validation folds (default = leaveone out)
%                       NOTE: if not using leave-one-out, computes average
%                       of 120 iterations (Varoquaux et al. 2017
%                       Neuroimage)
% r_method (optional)   : correlation method (1 = Pearson (default); 2 = spearman;
%                       3 = robust regress; 4 = Pearson partial using part_var;
%                       5 = Spearman partial using part_var          
% pthresh (optional)    : p threshold for feature selection (default = 0.01)
% part_var (optional)   : partial corr variable (leave blank if not using)
% motion_var (optional) : head motion as FD (if included, removes subjects with FD>0.15)
% outname (optional)    : name for output files (default = 'test')
% no_iter (optional)    : number of iterations for permute test (default=1000)
% OUTPUTS:
% p                     : p value for permutation test
% R_posneg              : real R value (mean of 100 iterations if using kfold)
% R_permute             : distribution of null R values

%% Settings
FD_thr=.20; % cutoff for removing subjects based on FD
if nargin<5 || isempty(r_method)
    r_method=1;
end
if nargin<7 || isempty(part_var)
    part_var=[];
end
if nargin<8 || isempty(motion_var)
    motion_var=[];
end

%% remove subjects with missing behavioral data
missing=isnan(all_behav);
all_mats(:,:,missing)=[];
if ~isempty(motion_var) 
    motion_var(missing)=[];
end
if ~isempty(part_var) 
    part_var(missing)=[];
end
all_behav(missing)=[];

%% remove high-motion subjects
if ~isempty(motion_var) 
    rm_subs=find(motion_var>FD_thr);
    if r_method~=3
        display(['removing ' num2str(length(rm_subs)) ' subjects due to high motion']);
        all_behav(rm_subs)=[];
        all_mats(:,:,rm_subs)=[];
    if ~isempty(part_var)
        part_var(rm_subs)=[];
    end
    end
    motion_var(rm_subs)=[];
end

%% Defaults
no_sub=length(all_behav);
if nargin<4 || isempty(kfolds)
    kfolds=no_sub;
end
if nargin<6 || isempty(pthresh)
    pthresh=0.01;
end
if nargin<9 || isempty(outname)
    outname='test';
end
if nargin<10 || isempty(no_iter)
    no_iter=1000;
end

if kfolds==no_sub % for leave-one-out, get one value
[R]=CPM_internal(all_mats,all_behav,dataset,kfolds,r_method,pthresh,part_var,motion_var,outname);
R_pos=R(1); R_neg=R(2); R_posneg=R(3);
else % for kfolds, do 120 iterations and get avg value
    for i=1:120
        [R]=CPM_internal(all_mats,all_behav,dataset,kfolds,r_method,pthresh,part_var,motion_var,outname); 
        R_pos(i)=R(1); R_neg(i)=R(2); R_posneg(i)=R(3);
    end
    R_pos=mean(R_pos); R_neg=mean(R_neg); R_posneg=mean(R_posneg);
end
% permutation test for significance for each subject
       R_permute_pos=[]; R_permute_neg=[]; R_permute_posneg=[];
        for it=2:no_iter
           display(['Performing iteration ' num2str(it)]); 
           % Permute labels
           if r_method==1 || r_method==2
            new_behav = all_behav(randperm(length(all_behav)));
            new_partial=[];
           elseif r_method==4 || r_method==5
            order = 1:length(all_behav);
            new_order = order(randperm(length(order)));
            new_behav = all_behav(new_order);
            new_partial = part_var(new_order);
           end
            [R_shuffled]=CPM_internal(all_mats,new_behav,'rand',kfolds,r_method,pthresh,part_var,motion_var);
            R_permute_pos=[R_permute_pos; R_shuffled(1)];
            R_permute_neg=[R_permute_neg; R_shuffled(2)];
            R_permute_posneg=[R_permute_posneg; R_shuffled(3)];
        end
        
% assess significance 
true_prediction_r_posneg=R_posneg;
prediction_r_posneg=[true_prediction_r_posneg; R_permute_posneg];
sorted_prediction_r_posneg=sort(prediction_r_posneg,'descend');
position_true_posneg=find(sorted_prediction_r_posneg==true_prediction_r_posneg);
p_posneg=position_true_posneg/no_iter;

true_prediction_r_pos=R_pos;
prediction_r_pos=[true_prediction_r_pos; R_permute_pos];
sorted_prediction_r_pos=sort(prediction_r_pos,'descend');
position_true_pos=find(sorted_prediction_r_pos==true_prediction_r_pos);
p_pos=position_true_pos/no_iter;

true_prediction_r_neg=R_neg;
prediction_r_neg=[true_prediction_r_neg; R_permute_neg];
sorted_prediction_r_neg=sort(prediction_r_neg,'descend');
position_true_neg=find(sorted_prediction_r_neg==true_prediction_r_neg);
p_neg=position_true_neg/no_iter;

p_all=[p_pos p_neg p_posneg];
R_permute_all=[R_permute_pos, R_permute_neg, R_permute_posneg];
