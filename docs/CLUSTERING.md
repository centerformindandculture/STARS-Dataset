# Contextual clustering and mortality comparisons

County clustering includes all 56 Montana counties and eight contextual features:
median age, American Indian/Alaska Native population percentage, median household
income, poverty, one-person households, divorced/separated percentage, associate
degree attainment, and uninsured percentage. Month clustering uses eight weather
and crime measures; weekday clustering uses crash, substance-mortality and happiness
measures. Suicide outcomes and decedent characteristics are excluded from clustering.

Features are standardized to z-scores. K-means uses 100 random starts. Diagnostics
include within-cluster sum of squares (elbow), mean silhouette and gap statistics
with 200 simulations. Selection is constrained to k=3–4: the smallest k with
silhouette within 0.02 of the best in that range. Unconstrained diagnostic optima
are reported separately and need not agree. Both k=3 and k=4 solutions are retained.
Labels are ordered by county income, month temperature or weekday happiness.
Leave-one-unit-out stability refits use the full-data feature scaling.

County-year Poisson models include cluster and categorical year. Count comparisons
have no population offset; rate comparisons include log population. Covariance
clusters the score contributions across observations within each county and is
multiplied by G/(G-1) * (n-1)/(n-p), with G=56 counties, n model rows and p parameters.
Coefficient intervals use t(55); omnibus comparisons use F(q,55). Counties are
assumed independent. Crude pooled rate intervals are exact Poisson intervals with
fixed denominators, separate from the year-adjusted model tests.

Age-standardized rates cover 2012–2021. Age-and-sex-standardized rates cover
2013–2021 because sex is unavailable in 2012. Models include cluster, demographic
stratum, categorical year and log population. Fitted rates are averaged using the
common pooled Montana stratum-by-year population distribution. Delta-method
uncertainty on the log-rate scale uses the county covariance and t(55) limits.
The standard is the study population, not a national standard population.

Monthly mortality comparisons use 2018, 2020 and 2021 (875 deaths); weekday
comparisons use 2013–2016, 2018, 2020 and 2021 (1,872 deaths). These years have
complete recording of the respective date variable. Offsets combine statewide
population and calendar-day exposure. Covariance is Poisson covariance multiplied
by max(1, Pearson dispersion), with residual-df t/F inference.

Holm adjustment covers four selected-solution omnibus tests: county counts, county
rates, month rates and weekday rates. The two demographic-standardization tests
form a separate family. Pairwise comparisons are adjusted within each fitted model.
Intervals are pointwise. Cluster-defining features receive no post-clustering tests.

Outputs are in `outputs/clustering/geographic`, `monthly`, `weekday` and `overview`.
The overview holds joint selection, stability, model and multiplicity summaries.
Population-estimate error, cluster-selection uncertainty and spatial dependence
between counties are not incorporated in the model intervals.
