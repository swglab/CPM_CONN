# CPM_CONN
Connectome-based predictive modeling analysis with CONN toolbox outputs

Setup:

1. In startup.m file, specify the parent directory of your project/dataset folders. Example (change to your specific directory):
   global globalDataDir;
   globalDataDir='/work/swglab/Aaron/data';
2. Preprocess your dataset in Conn, including the extraction of ROIs from an atlas (see utils) selected as an "atlas file" within Conn
3. Extract head motion (frame-wise displacement) in Conn (Setup > Covariates 1st level > Covariate tools > Compute new/derived first-level covariates > Compute 'FD_jenkinson')
   
   
Functions:

CPM_internal.m runs CPM within a dataset using an input vector for the behavioral variable and an input 3D matrix for functional connectivity (node x node x subject).

CPM_view_networks.m can be used to view intra- and inter-network contributions to a pre-computed CPM positive and negative masks
