function [p_allsubs,p_group_signrank]=CPM_internal_permute(all_mats,all_behav,subs,suffix,no_iterations,corr_method,partial_var)

% Inputs:
% 1. all_mats = ROI x ROI x trials FC matrix
% 2. all_behav = behavioral score vector
% 3. subs = list of subjects
% 4. suffix = suffix to label output figures
% 5. no_iteration = number of iterations for permutation test
% 6. corr_method: 1 = Pearson correlation (default), 2 = spearman, 4 =
% partial corr (Pearson) on lefout data 
% correlation, 3 = robust regression, 4 = partial corr
% 7. partial_var = vector of variabile to control for with partial corr
% method (e.g. FD)

%% Defaults
global globalDataDir;
datapath=[globalDataDir filesep 'MW_gradCPT'];
load('cdcol.mat');
cd([datapath]);
cd([datapath]);
if nargin<3 || isempty(subs)
    subs=importdata('sublist_HC.txt');
end
if nargin<5 || isempty(no_iterations)
   no_iterations=1000; 
end
if nargin<6 || isempty(corr_method)
   corr_method=1; 
end
if nargin<7 || isempty(partial_var)
   partial_var=[];
end

[R_posneg,R_pos,R_neg]=MW_CPM(all_mats,all_behav,subs,suffix,1,corr_method,partial_var);

% permutation test for significance for each subject
       R_permute=[]; 
        for it=2:no_iterations
           display(['Performing iteration ' num2str(it)]); 
           % Permute labels
           if corr_method==1 || corr_method==2
           new_behav = all_behav(randperm(length(all_behav)));
           new_partial=[];
           elseif corr_method==4 || corr_method==5
           order = 1:length(all_behav);
           new_order = order(randperm(length(order)));
           new_behav = all_behav(new_order);
           new_partial = partial_var(new_order);
           end
           [R_posneg_shuffled]=MW_CPM(all_mats,new_behav,subs,'rand',1,corr_method,new_partial);
           R_permute=[R_permute; R_posneg_shuffled];
        end
        
% assess significance for each subject
for i=1:size(R_permute,2)
    true_prediction_r(i)=R_posneg(i);
    prediction_r=[true_prediction_r(i); R_permute(:,i)];
    mean_null(i)=mean(R_permute(:,i));
    sorted_prediction_r=sort(prediction_r(:,1),'descend');
    position_true=find(sorted_prediction_r==true_prediction_r(i));
    p_allsubs(i)=position_true(1)/no_iterations;
end

% group-level t-test
[h,p_group_ttest]=ttest(true_prediction_r,mean_null);
p_group_signrank=signrank(true_prediction_r,mean_null);

% save p-value
cd([datapath filesep 'Group' filesep 'Results']);
save(['p_signrank_' suffix],'p_group_signrank');

% plot subjects
cd([datapath filesep 'Group' filesep 'Results']);
mkdir('figs'); cd('figs');
[r_plot,order]=sort(true_prediction_r,'descend');
p_plot=p_allsubs(order);
null_plot=mean_null(order);
%ru_plot=ru_posneg(order); rl_plot=rl_posneg(order);
color_plot=repmat(cdcol.white,length(r_plot),1);
color_null=repmat(cdcol.black,length(null_plot),1);

mean_r_allsubs=mean(true_prediction_r);
mean_null_allsubs=mean(mean_null);
ste_r=std(true_prediction_r)/sqrt(length(true_prediction_r));
ste_null=std(mean_null)/sqrt(length(mean_null));
color_plot=[color_plot; .7 .7 .7];
color_null=[color_null; cdcol.black];
r_plot=[r_plot mean_r_allsubs];
null_plot=[null_plot mean_null_allsubs];
% for i=1:length(r_plot)
%    if p_plot(i)<0.05
%        color_plot(i,:)=cdcol.grey;
%    end
% end

FigHandle=figure('Position',[200 200 700 400])
b=bar(r_plot,.4);
b.FaceColor='flat';
b.CData(:,:)=color_plot;
set(gca,'FontSize',13,'LineWidth',.5,'TickDir','in','box','off');
ylabel({['{\it r}-value:'];['Predicted vs. observed mind-wandering']});
hold on;
c=bar(null_plot,.4); c.FaceColor='flat'; c.CData(:,:)=color_null;
hold on;
errorbar(length(r_plot),r_plot(end),ste_r,'Color',cdcol.black);
hold on;
errorbar(length(null_plot),null_plot(end),ste_null,'Color',[.7 .7 .7]);
hold on;
%legend({['test'];['test2'];['test3']});
yline(0,'k--')
set(gca,'xtick',[])
xticks(length(r_plot));
xticklabels('Mean');
set(gca,'xcolor',[1 1 1])
ax=gca;
ax.XTickLabel{1}=['\color{black}' ax.XTickLabel{1}];
xlabel('Individual Participants','Color',cdcol.black);
if p_group_signrank<0.05
    sigstar_NoLine([length(r_plot),length(r_plot)])
end
%pause; close;
print(['predicted_MW_allsubs_internal_' suffix],'-r600','-dpng');
pause; close;

% FigHandle=figure('Position',[200 200 700 400])
% e=errorbar(1:length(r_plot),r_plot,(r_plot-rl_plot),(r_plot-ru_plot))
% e.Marker='.'; e.MarkerSize=10; e.Color=cdcol.black; e.CapSize=15;
% set(gca,'FontSize',12,'LineWidth',.5,'TickDir','in','box','off');
% title(['p=' num2str(p_group_signrank)],'Fontweight','normal');
% ylabel('Within-subject r-value');
% xlabel('Individual participants');
% yline(0,'k--')
% hold on;
% c=scatter(1:length(null_plot),null_plot,20,cdcol.scarlet,'filled');
% set(gca,'xtick',[])
% set(gca,'xcolor',[1 1 1])
% print(['predicted_MW_allsubs_internal_CIs_' suffix],'-r600','-dpng');
%pause; close;
