Guide to organization of this R project:

Getting the data:
- The folders 'raw_data/', 'intermediate_data/', and 'modeling_data/' are set to gitignore
- In order to reproduce results, you should have a copy of or link to Zongqi's Box folder within the 'raw_data' folder of the project and then running preprocess.R with the variable readin = FALSE. (Note: once this script is run once in its entirety, you can set readin = TRUE to bypass preliminary data cleaning and aggregation, which tends to take the most time).

Script overview
1. preprocess.R
2. tune_parameters.R - produces ROC plots for different values of tp and tw using Lasso cv.glmnet
-- (model.R) fits a Lasso model and examines the coefficients
-- (threshold.R) calculates relevant statistics for the chosen model
3. Other models (scripts/modeling/)
-- Adaptive lasso (adaptive_lasso.R)
-- Random forest (trees.R)
-- PCA (pca.R)
-- Sparse PCA (sparse_pca.R)
4. Plotting scripts (scripts/plotting/)
-- Probability of relapse over time per patient in test set (prob_plots.R)
-- Diagnostic Heatmaps of covariates over time for selected misclassified test patients (heatmap.R)
-- Correlation plots of predictors (corr_plot.R)


