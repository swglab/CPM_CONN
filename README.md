# CPM_CONN
Connectome-based predictive modeling analysis with CONN toolbox outputs

**Setup** (before using this code):

1. Preprocess your dataset in Conn, including extraction of ROIs from an atlas file (see utils) selected as an "atlas file" within Conn
2. Extract head motion as frame-wise displacement in Conn (Setup > Covariates 1st level > Covariate tools > Compute new/derived first-level covariates > Compute 'FD_jenkinson')
3. In startup.m file, specify the parent directory of your dataset folders. Example (change to your specific directory):
   global globalDataDir;
   globalDataDir='/work/swglab/Aaron/data';
4. Create a .mat cell array file with a list of subject names included in your Conn project
5. Create a .mat file with a vector of behavioral scores for each subject
   
**Functions:**

*extract_CONN_atlas_FC.m*: extracts functional connectivity matrices (from atlas) and mean FD, then merges across subjects (for input to CPM_internal.m)

*CPM_internal.m*: runs CPM within a dataset (kfold, leave one out, or use entire dataset to define and save model parameters)

*CPM_internal_permute.m*: runs permutation test to assess significance

*CPM_external.m*: test CPM (defined by CPM_internal.m output) in external data

*CPM_view_networks.m*: view intra- and inter-network contributions to positive and negative edges of a pre-computed CPM 

![test](https://github.com/swglab/CPM_CONN/blob/master/images/eg_feature_contributions.png)

*univariate_SchaeferYeo.m*: correlate behavior vs. connectivity between all intra- and inter-network pairs in the Schaefer atlas; apply FDR correction
