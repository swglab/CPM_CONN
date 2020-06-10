function [beta_lme]=univariate_FC_SchaeferYeo_trials(all_mats,all_behav,sub_code,suffix)

% compute mean within- and between-network Yeo7 network FC (based on
% Schaefer nodes), then perform linear mixed effects analysis to predict behavior
% written by Aaron Kucyi, Northeastern University

%% Defaults
global globalDataDir; % this must be set within startup.m
datapath=[globalDataDir];
net_names={'DMN'; 'DAN'; 'SAL'; 'FPCN'; 'LIM'; 'VIS'; 'SMN'};
global globalMaskDir;
net_labels=importdata([globalMaskDir filesep 'Schaefer' filesep 'net_label300.mat']);
mat_colors=cbrewer('div','RdBu',256);

%% compute within/between network FC and normalize within subjects
no_node=size(all_mats,1);
no_sub=max(sub_code);
no_trials_all=length(all_behav);
no_trials=length(all_behav)/no_sub;
sub_trials=1:size(all_mats,3);
trial_FC_norm_allsubs=[];
curr_sub_trials=1;
for i=1:max(sub_code)
    curr_trials = []; curr_mat = []; trial_FC=[];
    trial_FC_norm=NaN(length(net_names),length(net_names),no_trials);
    curr_trials = sub_trials(curr_sub_trials:curr_sub_trials+no_trials-1);
    curr_mat=all_mats(:,:,curr_trials);
    for j=1:max(net_labels) % seed network
        curr_network= [];
        curr_network=curr_mat(:,find(net_labels==j),:);
        for k=1:max(net_labels) % target network
           curr_target=[];
           curr_target=curr_network(find(net_labels==k),:,:);
           for l=1:size(curr_target,3) % loop through trial
               curr_trial=[];
               curr_trial=squeeze(curr_target(:,:,l));
               curr_trial=curr_trial(:);
               trial_FC(j,k,l)=mean(curr_trial(~isinf(curr_trial)));
           end
        end
    end
    % normalize FC within subject
    for j=1:size(trial_FC,1)
        for k=1:size(trial_FC,2)
            curr_edge=[];
            curr_edge=zscore(squeeze(trial_FC(j,k,:)));
            trial_FC_norm(j,k,:)=curr_edge;
        end
    end
    % concatenate across subjects
    trial_FC_norm_allsubs=cat(3,trial_FC_norm_allsubs,trial_FC_norm);
    curr_sub_trials=curr_sub_trials+no_trials;
end

% LME for each intra- and inter-network pair vs behavior
lmm_behav=all_behav;
   for i=1:size(trial_FC_norm_allsubs,1)
       for j=1:size(trial_FC_norm_allsubs,2)
           curr_edge=[];
           curr_edge=squeeze(trial_FC_norm_allsubs(i,j,:));
           corr_mat(i,j)=corr(all_behav,curr_edge);
           % LME
           lmm_variables=[]; lmm_table=[];
           lmm_variables=[sub_code all_behav curr_edge];
           lmm_table=table(lmm_variables(:,1),lmm_variables(:,2),lmm_variables(:,3),...
            'VariableNames',{'subject','behav','FC'});
            lme=fitlme(lmm_table,'behav~FC+(1|subject)');
            lme_anova=anova(lme,'DFmethod','satterthwaite');
            beta=lme.fixedEffects; beta=beta(2);
            p=lme_anova.pValue; p=p(2);
            beta_lme(i,j)=beta;
            p_lme(i,j)=p;
            F=lme_anova.FStat; F=F(2);
            F_lme(i,j)=F;
       end
   end
   p_lme=tril(p_lme);
   p_lme(p_lme==0)=NaN;
   beta_lme=tril(beta_lme);
   %beta_lme(beta_lme==0)=NaN;

%% fdr correction
p=p_lme(:); p(isnan(p))=[]; p_fdr=fdr(p);
p_mat=isnan(p_lme); p_fdr_lme=NaN(size(p_lme,1),size(p_lme,2));
for i=1:length(p_fdr)
    ind=[];
    ind=find(p_mat==0);
    p_fdr_lme(ind(1))=p_fdr(i);
    p_mat(ind(1))=1;
end
[y_fdr,x_fdr]=find(p_fdr_lme<0.05); % coordinates for plotting *

%% plot beta weights
mat_colors(128,:)=[1 1 1];
mat_colors(129,:)=[1 1 1];
cd([savepath]); mkdir('figs'); cd('figs');
anchor=max(beta_lme(:))+.05;

fig=figure('Position',[100 100 600 600]);
nx=length(net_names); %ny=length(net_names_plot);
imagesc(beta_lme,[-anchor anchor]);
colormap(flipud(mat_colors)); 
c=colorbar('Location','southoutside');
c.FontSize=20;
c.Label.String={'\beta'};
hold on;
cut=1; cut2=-1;
for i=1:nx+1
    plot([0.5,cut+0.5],[i-.5,i-.5],'k-');
    plot([i-.5,i-.5],[0.5+cut2,0.5+nx],'k-');
    cut=cut+1; cut2=cut2+1;
end
set(gca,'xtick',1:size(beta_lme,2))
box off
set(gca,'xticklabel',[net_names(1:length(net_names))])
set(gca,'ytick',1:size(beta_lme,2))
set(gca,'yticklabel',[net_names(1:length(net_names))])
set(gca,'TickLength',[0 0])
set(gca,'FontSize',20,'FontName','Arial');
xtickangle(45);
hold on;
s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor','k')
legend(s,['{\it p_{FDR}} <0.05']);
legend boxoff
title({['Univariate Analysis'];[' ']},'Fontweight','normal');
print(['Univariate_SchaeferYeo_vs_MW_' suffix],'-r600','-dpng');  


pause; close;