# two_phase_AUC
R code implementing estimators of Area Under the Curve(AUC) in two-phase studies

The "AUC_code.R" file contains the main R function "compute_AUC()" that implements three estimators of AUC and reports the corresponding point estimates and the variances: 

(i) The full cohort estimator is based on the partial risk factor information coming from all the subjects in the cohort     (Phase I subjects). The Delong estimator of variance is implemented.

(ii) The Inverse Probability Weighted (IPW) estimator is based on the second phase subjects who contribute full risk factor information. The influence function based variance estimator (reported in the paper) is implemented.

(iii) The Two-Phase Estimator is the estimator proposed in the paper. It uses information from the second phase subjects with complete risk factor data and the subjects not included in the second phase with partial risk factor data. The influence function based variance estimator (reported in paper) is implemented.

The compute_AUC() function takes the following arguments:

cohort.data: The dataset from Phase-I subjects with partial risk factor information. For illustrative purposes, I include one simulated dataset with columns: 
                 id - Subject identifier,
                 X1 - Continuous risk-factor,
                 X2 - Continuous risk-factor,
                 X3 - Continuous risk-factor,
                 X4 - Continuous risk-factor,
                 X5 - Binary risk-factor,
                 X6 - Binary risk-factor,
                 study_entry_age - Age (in years) of entry of a subject to the study,
                 study_followup - Length of followup (in years) of a subject,
                 study_exit_age - Age (in years) of exit of a subject from the study,
                 time_to_event - Time (in years) of onset of the disease, for censored subjects it is Inf,
                 event - binary indicator of disease coded as 1 for event and 0 for censored,
                 obs.riskscore - risk-score based on the risk factors X1,...,X6,
                 sampling.weights - selection probability of a Phase-I subject to Phase-II,
                 include - binary indicator of inclusion of a Phase-I subject to Phase-II
              
case.control.data: The dataset from Phase-II subjects (subset of Phase I subjects included in Phase II) with complete risk factor information. For illustrative purposes, I include one simulated dataset. In addition to the columns described above,it includes two additional binary risk-factors X7 and X8. Only the Phase-II subjects have information on X7 and X8.

full.model.formula: Model formula for the specification of the full risk-score (i.e., linear predictor associated with the full set of risk factors). The Phase II subjects have information on the full set of risk-factors

reduced.model.formula: Model formula for the specification of the partial risk-score (i.e., linear predictor associated with the partial set of risk factors). The Phase-I subjects not include in Phase II have information on the partial set of risk factors

risk.factors.logRR: Log-relative risks associated with the risk-factors in the full model, known apriori from earlier studies and used to construct the risk scores

no.riskscore.cat: Number of categories used for binning the partial risk-scores (e.g., for deciles, set the value at 10)

The compute_AUC() function returns the following objects:

Full_Cohort_Estimator: The full cohort estimator is based on the partial risk factor information coming from all the subjects in the cohort (Phase I subjects). 

Variance_Full_Cohort_Estimator: Delong variance formula for the variance of the full cohort estimator

IPW_Estimator: The Inverse Probability Weighted (IPW) estimator is based on the second phase subjects who contribute full risk factor information

Variance_IPW_Estimator: The influence function based variance estimator of the IPW estimator

Two_Phase_Estimator: The Two-Phase Estimator is the estimator proposed in the paper. It uses information from the second phase subjects with complete risk factor data and the subjects not included in the second phase with partial risk factor data.

Variance_Two_Phase_Estimator: The influence function based variance estimator of the two-phase estimator




  
