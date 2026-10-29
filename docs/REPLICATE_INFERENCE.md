# Regression uncertainty

Point estimates use mortality-record weight 1 and ACS person weight PERWT.
Annual person weights are not divided by the number of study years. The 80
REPWTP1–REPWTP80 weights estimate ACS sampling covariance; REPWTP is a presence flag.

For each supported year y, each replicate r replaces PERWT only for that year's
ACS records. Other annual weights and all mortality weights remain fixed. With
full coefficient vector b and replicate vector b_yr:

    V_ACS = sum_y (4/80) sum_r (b_yr-b)(b_yr-b)'

Annual ACS samples are assumed independent. This procedure does not estimate
cross-year sampling covariance. Signed replicate weights are retained. Replicate
fits solve the logistic score equations directly, aggregating identical categorical
design rows into covariate/year cells. For case counts D, population totals C,
design matrix X and p=logistic(Xb), the score is X'[D-(D+C)p]=0.

Let H=X' diag((D+C)p(1-p))X. Mortality-event covariance is

    V_case = H^-1 [X' diag(D(1-p)^2)X] H^-1

The reported covariance is V_case+V_ACS, assuming independence between ACS sampling
error and mortality-event variation. Intervals use normal Wald limits on the
log-odds scale. Coefficient tests, significance stars, within-sex age contrasts and
joint age-by-sex Wald tests use this covariance. The model-based covariance is
reported as a methodological comparator, not added to the two components.

Checks cover score residuals, finite estimates, information matrices, full-sample
reconstruction and agreement with binomial fitting under positive-weight
perturbation. Regression diagnostics contain `replicate_model_checks.csv`,
`replicate_weight_summary.csv`, and `covariance_method_comparison.csv`.
Each analysis stores aggregate replicate coefficient vectors and the case-process,
ACS and combined covariance matrices under `replicate_details`.

These intervals do not account for predictor classification errors, missingness,
unmeasured confounding or dependence among deaths. See [IPUMS replicate-weight
guidance](https://usa.ipums.org/usa/repwt.shtml) and [Census ACS PUMS documentation](https://www.census.gov/programs-surveys/acs/microdata/documentation.html).
