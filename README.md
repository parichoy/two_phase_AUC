# two_phase_AUC
R code implementing estimators of Area Under the Curve(AUC) in two-phase studies

For running the example code, one needs to run the R script file "AUC_code.R". The code is designed to run in a way that the three files: "cohort.data1.txt", "case.control.data1.txt" and "AUC_code.R" are in the same local directory.

The "AUC_code.R" file contains the main R function "compute_AUC()" that implements three estimators of AUC and reports the corresponding point estimates and the variances: 

(i) The full cohort estimator is based on the partial risk factor information coming from all the subjects in the cohort (Phase I subjects). The Delong estimator of variance is implemented.

(ii) The Inverse Probability Weighted (IPW) estimator is based on the second phase subjects who contribute full risk factor information. The influence function based variance estimator (reported in the paper) is implemented.

(iii) The Two-Phase Estimator is the estimator proposed in the paper. It uses information from the second phase subjects with complete risk factor data and the subjects not included in the second phase with partial risk factor data. The influence function based variance estimator (reported in paper) is implemented.

The compute_AUC() function takes the following arguments:

cohort.data1: The dataset from Phase-I subjects with partial risk factor information. For illustrative purposes, I include one simulated dataset with columns: 
                 id - Subject identifier,
                 X1 - Continuous risk-factor,
                 X2 - Continuous risk-factor,
                 X3 - Continuous risk-factor,
                 X4 - Continuous risk-factor,
                 X5 - Binary risk-factor,
                 X6 - Binary risk-factor,
                 study_entry_age - Age (in years) of entry of a subject to the study,
                 followup - Length of followup (in years) of a subject,
                 study_exit_age - Age (in years) of exit of a subject from the study,
                 time_to_event - Time (in years) of onset of the disease, for censored subjects it is Inf,
                 event - binary indicator of disease coded as 1 for event and 0 for censored,
                 sampling.weights - selection probability of a Phase-I subject to Phase-II,
                 include - binary indicator of inclusion of a Phase-I subject to Phase-II
              
case.control.data1: The dataset from Phase-II subjects (subset of Phase I subjects included in Phase II) with complete risk factor information. For illustrative purposes, I include one simulated dataset. In addition to the columns described above, it includes two additional binary risk-factors X7 and X8. Only the Phase-II subjects have information on X7 and X8.

full.model.formula: Model formula for the specification of the full risk-score (i.e., linear predictor associated with the full set of risk factors). The Phase II subjects have information on the full set of risk-factors

reduced.model.formula: Model formula for the specification of the partial risk-score (i.e., linear predictor associated with the partial set of risk factors). The Phase-I subjects not included in Phase II have information on the partial set of risk factors

risk.factors.logRR: Log-relative risks associated with the risk-factors in the full model, known apriori from earlier studies and used to construct the risk scores

no.riskscore.cat: Number of categories used for binning the partial risk-scores (e.g., for deciles, set the value at 10)

The compute_AUC() function returns the following objects:

Full_Cohort_Estimator: The full cohort estimator is based on the partial risk factor information coming from all the subjects in the cohort (Phase I subjects). 

Variance_Full_Cohort_Estimator: Delong variance formula for the variance of the full cohort estimator

IPW_Estimator: The Inverse Probability Weighted (IPW) estimator is based on the second phase subjects who contribute full risk factor information

Variance_IPW_Estimator: The influence function based variance estimator of the IPW estimator

Two_Phase_Estimator: The Two-Phase Estimator is the estimator proposed in the paper. It uses information from the second phase subjects with complete risk factor data and the subjects not included in the second phase with partial risk factor data.

Variance_Two_Phase_Estimator: The influence function based variance estimator of the two-phase estimator



Estimators of AUC accounting for differential follow-up of subjects

The "AUC_code.R" also includes a function "compute_AUC_adj" that implements estimators of AUC for a model predicting risk over a fixed time interval (e.g., 10 years) accounting for differential follow-up of subjects.

We create two simulated datasets "cohort.data2.txt" and "case.control.data2.txt". The first dataset includes Phase-I subjects who have either developed disease by 10-years or were followed up for at least 10 years. The second dataset includes only those Phase-I subjects who developed the disease by 10-years or were followed up for at least 10 years and included in the second phase.

For running this part of the code, one needs to run the R script file "AUC_code.R". The code is designed to run in a way that the three files: "cohort.data2.txt", "case.control.data2.txt" and "AUC_code.R" are in the same local directory.

The "AUC_code.R" file contains the main R function "compute_AUC_adj()" that implements three estimators of AUC and reports the corresponding point estimates and the variances: 

(i) The adjusted full cohort estimator is based on the partial risk factor information coming from the subjects in the cohort who have developed the disease by 10-years or were followed up for at least 10-years (Phase I subjects). A variance estimator is also implemented that accounts for the additional uncertainty due to selection of subjects by using the inverse of the probability that a control subject is followed up for at least 10-years.

(ii) The Inverse Probability Weighted (IPW) estimator is based on the second phase subjects who contribute full risk factor information. The influence function based variance estimator (reported in the paper) is implemented that accounts for both the selection of subjects in the Phase-I dataset and selection of these Phase-I subjects to the Phase-II sample using an inverse probability weighted approach.

(iii) The Two-Phase Estimator is the estimator proposed in the paper. It uses information from the second phase subjects with complete risk factor data and the subjects not included in the second phase with partial risk factor data. The influence function based variance estimator (reported in paper) is implemented that also accounts for the additional due to selection of  Phase-I subjects by using the inverse of the probability that a control subject is followed up for at least 10-years..

The compute_AUC_adj() function takes the following arguments:

cohort.data2: The dataset from Phase-I subjects who have developed the disease by 10-years or were followed up for at least 10-years with partial risk factor information. For illustrative purposes, I include one simulated dataset with columns: 
                 id - Subject identifier,
                 X1 - Continuous risk-factor,
                 X2 - Continuous risk-factor,
                 X3 - Continuous risk-factor,
                 X4 - Continuous risk-factor,
                 X5 - Binary risk-factor,
                 X6 - Binary risk-factor,
                 study_entry_age - Age (in years) of entry of a subject to the study,
                 followup - Length of followup (in years) of a subject,
                 sampling.weights.phase1 - selection probability of a subject (i.e., a case within 10-years or follow-up at least 10-years) to the Phase-I sample. 
                 study_exit_age - Age (in years) of exit of a subject from the study,
                 time_to_event - Time (in years) of onset of the disease, for censored subjects it is Inf,
                 event - binary indicator of disease coded as 1 for event and 0 for censored,
                 sampling.weights - selection probability of a Phase-I subject to Phase-II,
                 include - binary indicator of inclusion of a Phase-I subject to Phase-II
              
case.control.data: The dataset from Phase-II subjects (subset of Phase I subjects included in Phase II) with complete risk factor information. For illustrative purposes, I include one simulated dataset. In addition to the columns described above,it includes two additional binary risk-factors X7 and X8. Only the Phase-II subjects have information on X7 and X8.

full.model.formula: Model formula for the specification of the full risk-score (i.e., linear predictor associated with the full set of risk factors). The Phase II subjects have information on the full set of risk-factors

reduced.model.formula: Model formula for the specification of the partial risk-score (i.e., linear predictor associated with the partial set of risk factors). The Phase-I subjects not included in Phase II have information on the partial set of risk factors

risk.factors.logRR: Log-relative risks associated with the risk-factors in the full model, known apriori from earlier studies and used to construct the risk scores

no.riskscore.cat: Number of categories used for binning the partial risk-scores (e.g., for deciles, set the value at 10)

The compute_AUC_adj() function returns the following objects:

Adjusted_Full_Cohort_Estimator: The full cohort estimator is based on the partial risk factor information coming from all the subjects in the cohort (Phase I subjects). 

Variance_Adjusted_Full_Cohort_Estimator: Delong variance formula for the variance of the full cohort estimator

Adjusted_IPW_Estimator: The Inverse Probability Weighted (IPW) estimator is based on the second phase subjects who contribute full risk factor information

Variance_Adjusted_IPW_Estimator: The influence function based variance estimator of the IPW estimator

Adjusted_Two_Phase_Estimator: The Two-Phase Estimator is the estimator proposed in the paper. It uses information from the second phase subjects with complete risk factor data and the subjects not included in the second phase with partial risk factor data.

Variance_Adjusted_Two_Phase_Estimator: The influence function based variance estimator of the two-phase estimator
