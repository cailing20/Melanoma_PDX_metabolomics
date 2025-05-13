incoming_data_QC.R:
Examine data structure, signal correlation with run order in QC samples, missingness of values and their patterns, record metabolites that show association between missingness and run order (non-random missingness), as ro_metab.rds

data_cleaning.R:
Using samples within interquartile range to perform correlation with run order, this identifies metabolites with signal significantly correlated with run order. For these metabolites, we then perform LOESS correction based on IQR samples. A set of sstable metabolites with signal between 10^5 to 10^10, no missing values and RSD < 5% were identified and used to construct a normalizing factor (surrogate for total ion count normalization). Data was log10 transformed and saved as work/cleaned_dat.rds, sample information saved as work/sample_info.rds

exploratory_analyses.R
Corrected the sample mislabeling problem, removed unreliable metabolites. Created imputed data
