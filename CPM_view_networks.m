function CPM_view_networks(cpm,atlas,dataset,norm)

% Extract top intra- and inter-network contributors to a CPM
% written by Aaron Kucyi, Northeastern University
% INPUTS:
% cpm (required)        : structure containing pos_mask and neg_mask
% atlas (required)      : 1=Shen268, 2=Schaefer300 (7Networks), 3=Schaefer300 (17Networks)
% dataset (required)    : name of dataset folder name
% norm (optional)       : 0 = no normalization; 1 = normalize # of features per network pair by total # of network pairs
% OUTPUTS:
% saves figure within cpm_results/figs

%% Settings
nfeatures=5; % # of top feature types to extract
pos_colors=cbrewer('seq','Reds',256);
neg_colors=cbrewer('seq','Blues',256);

global globalDataDir;
datapath=[globalDataDir];
if nargin<5 || isempty(norm)
   	norm=0;
end
%% Load network labels
global globalMaskDir
if atlas==1
atlas_name='Shen268';
net_labels=importdata([globalMaskDir filesep 'Shen268' filesep 'Lookup_shen268']);
net_labels=net_labels(:,7);
net_names={'SMN';'[]'; 'CO'; 'AUD'; 'DMN';'[]'; 'VIS'; 'FPN';...
    'SAL'; 'SUB'; 'VAN'; 'DAN'; 'Unknown'};
elseif atlas==2
net_labels=importdata([globalMaskDir filesep 'Schaefer' filesep 'net_label300.mat']);
net_names={'DMN'; 'DAN'; 'SAL'; 'FPCN'; 'LIM'; 'VIS'; 'SMN'};
atlas_name='Schaefer300';
elseif atlas==3 
net_labels=importdata([globalMaskDir filesep 'Schaefer' filesep 'net_label300_17networks.mat']);
net_names={'VIS_{A}'; 'VIS_{B}'; 'SMN_{A}'; 'SMN_{B}'; 'DAN_{A}'; 'DAN_{B}'; 'SAL_{A}'; 'SAL_{B}';...
 'LIM_{A}'; 'LIM_{B}'; 'FPCN_{B}'; 'FPCN_{A}'; 'FPCN_{C}'; 'DMN_{A}'; 'DMN_{B}'; 'DMN_{C}'; 'TP'};
atlas_name='Schaefer17Networks';
end

%% get n nodes per networks and n network-network pairs
for i=1:length(net_names)
    node_count_net(i)=length(find(net_labels==i));
end
pair_count=NaN(length(net_names),length(net_names));
for i=1:length(net_names)
   for j=1:length(net_names)
       pair_count(i,j)=node_count_net(i)*node_count_net(j);
    end
end

%% Pos mask
[x_pos,y_pos]=find(cpm.pos_mask==1);
xy_pos_labels=[net_labels(x_pos) net_labels(y_pos)];

% count # of network-pair feature types
pos_count=NaN(length(net_names),length(net_names));
for i=1:length(net_names)
    for j=1:length(net_names)
        pos_count(i,j)=length(find(xy_pos_labels(:,1)==i & xy_pos_labels(:,2)==j));
    end
end
pos_count=tril(pos_count);

%% neg mask
[x_neg,y_neg]=find(cpm.neg_mask==1);
xy_neg_labels=[net_labels(x_neg) net_labels(y_neg)];

% count # of network-pair feature types
neg_count=NaN(length(net_names),length(net_names));
for i=1:length(net_names)
    for j=1:length(net_names)
        neg_count(i,j)=length(find(xy_neg_labels(:,1)==i & xy_neg_labels(:,2)==j));
    end
end
neg_count=tril(neg_count);

%% normalize by total number of network-network pairs
if atlas==1 % remove "unknown" networks for Shen268 atlas
pair_count_known=pair_count(1:end-1,1:end-1);
pos_count_known=pos_count(1:end-1,1:end-1); % delete counts for unknown networks
neg_count_known=neg_count(1:end-1,1:end-1); % delete counts for unknown networks
elseif atlas==2 || atlas==3
pair_count_known=pair_count;
pos_count_known=pos_count; % delete counts for unknown networks
neg_count_known=neg_count;    
end
if norm==1
    pos_count_norm=pos_count_known./pair_count_known;
    neg_count_norm=neg_count_known./pair_count_known;
else
    pos_count_norm=pos_count_known;
    neg_count_norm=neg_count_known;
end

%% top feature types: pos mask
xpos_top=[]; ypos_top=[];
temp_pos_count=pos_count_norm;
for i=1:nfeatures
    [val_pos_top(i),ind_pos_top(i)]=max(temp_pos_count,[],'all','linear');
    curr_ind=find(temp_pos_count==val_pos_top(i));
    [x,y]=find(temp_pos_count==temp_pos_count(curr_ind(1)));
    xpos_top(i)=x(1);
    ypos_top(i)=y(1);
    temp_pos_count(ind_pos_top(i))=NaN;
end

%% top feature types: neg mask
xneg_top=[]; yneg_top=[];
temp_neg_count=neg_count_norm;
for i=1:nfeatures
    [val_neg_top(i),ind_neg_top(i)]=max(temp_neg_count,[],'all','linear');
    curr_ind=find(temp_neg_count==val_neg_top(i));
    [x,y]=find(temp_neg_count==temp_neg_count(curr_ind(1)));
    xneg_top(i)=x(1);
    yneg_top(i)=y(1);
    temp_neg_count(ind_neg_top(i))=NaN;
end

%% plot matrix of relative feature contributions
% positive mask
plot_pos_count=pos_count_norm;
if atlas==1
plot_pos_count([2 6],:)=[]; plot_pos_count(:,[2 6])=[];
end
net_names_plot=net_names;
if atlas==1
net_names_plot([2 6 13])=[];
end
pos_colors=[1 1 1; pos_colors(1:end-1,:)];
if atlas<3
fig=figure('Position',[0 0 900 450]);
elseif atlas==3
fig=figure('Position',[0 0 1300 600]);
end
ax(1)=subplot(1,2,1);
nx=length(net_names_plot); 
imagesc(plot_pos_count);
colormap(ax(1),pos_colors); 
c=colorbar('Location','southoutside');
c.FontSize=14;
if norm==1
 c.Label.String='Percentage of Edges';   
elseif norm==0
c.Label.String='Number of Edges';
end
hold on;
cut=1; cut2=-1;
for i=1:nx+1
    plot([0.5,cut+0.5],[i-.5,i-.5],'k-');
    plot([i-.5,i-.5],[0.5+cut2,0.5+nx],'k-');
    cut=cut+1; cut2=cut2+1;
end
set(gca,'xtick',1:size(pos_count_norm,2))
box off
set(gca,'xticklabel',[net_names_plot(1:length(net_names_plot))])
set(gca,'ytick',1:size(pos_count_norm,2))
set(gca,'yticklabel',[net_names_plot(1:length(net_names_plot))])
set(gca,'TickLength',[0 0])
set(gca,'FontSize',12,'FontName','Arial');
xtickangle(90);
hold on;

% negative mask
ax(2)=subplot(1,2,2);
plot_neg_count=neg_count_norm;
if atlas==1
plot_neg_count([2 6],:)=[]; plot_neg_count(:,[2 6])=[];
end
net_names_plot=net_names;
if atlas==1
net_names_plot([2 6 13])=[];
end
neg_colors=[1 1 1; neg_colors(1:end-1,:)];

nx=length(net_names_plot); %ny=length(net_names_plot);
imagesc(plot_neg_count);
colormap(ax(2),neg_colors);
c=colorbar('Location','southoutside');
c.FontSize=14;
if norm==1
 c.Label.String='Percentage of Edges';   
elseif norm==0
c.Label.String='Number of Edges';
end
hold on;
cut=1; cut2=-1;
for i=1:nx+1
    plot([0.5,cut+0.5],[i-.5,i-.5],'k-');
    plot([i-.5,i-.5],[0.5+cut2,nx+0.5],'k-');
    cut=cut+1; cut2=cut2+1;
end
set(gca,'xtick',1:size(neg_count_norm,2))
box off
set(gca,'xticklabel',[net_names_plot(1:length(net_names_plot))])
set(gca,'ytick',1:size(neg_count_norm,2))
set(gca,'yticklabel',[net_names_plot(1:length(net_names_plot))])
set(gca,'TickLength',[0 0])
set(gca,'FontSize',12,'FontName','Arial');
xtickangle(90);
mkdir([datapath filesep dataset filesep 'cpm_results' filesep 'figs']);
savepath=[datapath filesep dataset filesep 'cpm_results' filesep 'figs'];
print([savepath filesep 'feature_contributions_' atlas_name],'-r600','-dpng');  
pause; close;


