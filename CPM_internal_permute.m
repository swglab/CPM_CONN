function [p,R_posneg,R_permute]=CPM_internal_permute(all_mats,all_behav,dataset,...
    kfolds,r_method,pthresh,part_var,motion_var,outname,no_iter)

% Permutation test for CPM_internal.m
% written by Aaron Kucyi, Northeastern University
% adapted from Shen et al. (2017 Nature Protocols)
% INPUTS:
% all_mats (required)   : ROI x ROI x trials FC matrix (or single vector for one ROI/edge)
% all_behav (required)  : behavioral score vector
% dataset (required)    : name of dataset folder name
% kfolds (optional)     : number of cross-validation folds (default = leave one out) 
% r_method (optional)   : correlation method (1 = Pearson (default); 2 = spearman;
%                       3 = robust regress; 4 = Pearson partial using part_var;
%                       5 = Spearman partial using part_var          
% pthresh (optional)    : p threshold for feature selection (default = 0.01)
% part_var (optional)   : partial corr variable (leave blank if not using)
% motion_var (optional) : head motion as FD (if included, removes subjects with FD>0.15)
% outname (optional)    : name for output files (default = 'test')
% no_iter (optional)    : number of iterations for permute test (default=1000)
% OUTPUTS:
% p

%% Settings
FD_thr=.15; % cutoff for removing subjects based on FD
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

[R_posneg]=CPM_internal(all_mats,all_behav,dataset,kfolds,r_method,pthresh,part_var,motion_var,outname);

% permutation test for significance for each subject
       R_permute=[]; 
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
            new_partial = partial_var(new_order);
           end
            [R_posneg_shuffled]=CPM_internal(all_mats,new_behav,'rand',kfolds,r_method,pthresh,part_var,motion_var);
            R_permute=[R_permute; R_posneg_shuffled];
        end
        
% assess significance 
true_prediction_r=R_posneg;
prediction_r=[true_prediction_r; R_permute];
sorted_prediction_r=sort(prediction_r('descend'));
position_true=find(sorted_prediction_r==true_prediction_r);
p=position_true/no_iter;

% for i=1:size(R_permute,2)
%     true_prediction_r(i)=R_posneg(i);
%     prediction_r=[true_prediction_r(i); R_permute(:,i)];
%     mean_null(i)=mean(R_permute(:,i));
%     sorted_prediction_r=sort(prediction_r(:,1),'descend');
%     position_true=find(sorted_prediction_r==true_prediction_r(i));
%     p_allsubs(i)=position_true(1)/no_iterations;
% end

% group-level t-test
% [h,p_group_ttest]=ttest(true_prediction_r,mean_null);
% p_group_signrank=signrank(true_prediction_r,mean_null);





