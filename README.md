# CPM_CONN
Connectome-based predictive modeling analysis with CONN toolbox outputs

Setup:

1. In your startup.m file, specify the parent directory of project/dataset folders. Example (change second line to your directory):
   global globalDataDir
   globalDataDir='/work/swglab/Aaron/data';
   
Functions:

CPM_internal.m runs CPM within a dataset using an input vector for the behavioral variable and an input 3D matrix for functional connectivity (node x node x subject).

CPM_view_networks.m can be used to view intra- and inter-network contributions to a pre-computed CPM positive and negative masks
