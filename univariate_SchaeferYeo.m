function [beta_lme]=univariate_SchaeferYeo(all_mats,all_behav,suffix)

% compute mean within- and between-network Yeo7 network FC (based on
% Schaefer nodes), then perform linear mixed effects analysis to predict behavior
% written by Aaron Kucyi, Northeastern University

%% Defaults
global globalDataDir; % this must be set within startup.m
datapath=[globalDataDir];
net_names={'DMN'; 'DAN'; 'SAL'; 'FPCN'; 'LIM'; 'VIS'; 'SMN'};
net_labels=importdata(['Schaefer_label300_7networks.mat']);
mat_colors=cbrewer('div','RdBu',256);

%% compute within/between network FC and normalize within subjects
no_node=size(all_mats,1);
no_sub=length(all_behav);
sub_FC_all=[];
for i=1:no_sub
    curr_mat = []; 
    sub_FC=NaN(length(net_names),length(net_names));
    curr_mat=all_mats(:,:,i);
    for j=1:max(net_labels) % seed network
        curr_network= [];
        curr_network=curr_mat(:,find(net_labels==j),:);
        for k=1:max(net_labels) % target network
           curr_target=[];
           curr_target=curr_network(find(net_labels==k),:,:);
           curr_target=curr_target(:);
           sub_FC(j,k)=mean(curr_target(~isinf(curr_target)));
        end
    end
    % concatenate across subjects
    sub_FC_all=cat(3,sub_FC_all,sub_FC);
end

% EDIT HERE
% LME for each intra- and inter-network pair vs behavior
lmm_behav=all_behav;
   for i=1:size(sub_FC,1)    
        for j=1:size(sub_FC,2)
            curr_edge=[];
            curr_edge=squeeze(sub_FC_all(i,j,:));
            corr_mat(i,j)=corr(all_behav,curr_edge);
            % LME
            lmm_variables=[]; lmm_table=[];
            lmm_variables=[all_behav curr_edge];
            lmm_table=table(lmm_variables(:,1),lmm_variables(:,2),lmm_variables(:,3),...
             'VariableNames',{'subject','behav','FC'});
%             lme=fitlme(lmm_table,'behav~FC+(1|subject)');
%             lme_anova=anova(lme,'DFmethod','satterthwaite');
%             beta=lme.fixedEffects; beta=beta(2);
%             p=lme_anova.pValue; p=p(2);
%             beta_lme(i,j)=beta;
%             p_lme(i,j)=p;
%             F=lme_anova.FStat; F=F(2);
%             F_lme(i,j)=F;
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