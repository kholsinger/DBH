Associated with the manuscript "Fusing tree-ring and forest inventory
data to infer influences on tree growth" (Evans et al. 2017. Ecosphere
________), we publish R code to implement our statistical model in the
public GitHub repository located at https://github.com/kholsinger/DBH.
We point readers to a Readme.md file found in this repository, which
contains the same information as follows here.

We describe four models reported in the manuscript: coupled
vs. uncoupled, and 24 monthly climate variables vs. 4 seasonal climate
variables as predictors (see Methods, Evans et al. 2017. Ecosphere
________). The workflow followed in our analysis is controlled by one
of two R driver files, one each for the model with 24 monthly climate
variables (gi-plus-dbh.R) vs. 4 seasonal climate variables
(gi-plus-dbhSeas.R). These R driver files accomplish the following
steps of the workflow

1. bring in data,
2. prepare those data in data
objects for analysis,
3. set the stage for Stan model analysis,
4. fit a Stan model, and
5. save Stan model output.

The R driver files contain, at the top, several controls (TRUE/FALSE)
on what will be done, including debugging, implementation of the
coupled vs. uncoupled model, model comparison (i.e., between coupled
vs. uncoupled model output), writing results to a text file, and
saving Stan MCMC output.  Regardless of these settings, execution of
the R driver will source two data preparation scripts:

1. "prepare-data.R", which reads in tree-ring data (all .rwl files
located in the folder plot-data) and climate covariate data
(PRISM_MCN.txt and vpd_MCC_cru1901_1913.txt) and
2. "dbh-process-data.R", which reads in diameter data (MCNallDBH4.csv)
and plot covariate data (MCNplotCovariatesNEW.csv). Note that the CRU
VPD data actually includes the time series 1901-2013.

If coupled is TRUE and uncoupled is FALSE, the R driver file
gi-plus-dbh.R will use "gi-plus-dbh.stan" for the Stan model file; if
coupled is FALSE and uncoupled is TRUE, the R driver file
gi-plus-dbh.R will use "gi-plus-dbh-uncoupled.stan" for the Stan model
file. Similarly, if coupled is TRUE and uncoupled is FALSE, the R
driver file gi-plus-dbhSeas.R (testing 4 seasonal climate predictors)
will use "gi-plus-dbhM2.stan" for the Stan model file; if coupled is
FALSE and uncoupled is TRUE, the R driver file gi-plus-dbhSeas.R will
use "gi-plus-dbh-uncoupledM2.stan" for the Stan model file. Parallel
models that test the effect of fire + thinning are found in the branch
"Fire". Readers unfamiliar with Stan are encouraged to consult the
Stan website http://mc-stan.org/ to get oriented to the basic
structure and syntax of modeling in Stan.  If write.results.file is
set to TRUE in the R driver file, summary statistics of parameters
monitored during Stan's MCMC simulation are appended to the model
output file "results-gi-plus-dbh.txt". If save is TRUE, samples are
saved from the MCMC chains in an Rsave file (named either
"results-gi-pls-dbh.Rsave" or "results-gi-pls-dbh-uncoupled.Rsave"),
according to MCMC settings specified at the top of the R driver
file. If compare is TRUE, a comparison is made the MCMC samples in the
coupled vs. uncoupled Rsave files (see the description of
leave-one-out [LOO] cross-validation, under Model comparison in the
Methods, Evans et al. 2017 Ecosphere __________). Results from that
comparison are written to the file "model-comparison.txt".

Additional models explored are found in the repository, including
gi-plus-dbhM2a.stan, gi-plus-dbhM2b.stan, gi-plus-dbhM2c.stan,
gi-plus-dbhProfCurv.stan (testing the effect of the GIS-derived
topographic variable profile curvature), and gi-plus-dbhTRMI.stan
(testing the effect of a topographic relative moisture index at a
variety of spatial scales; citation). Even more models considered but
not reported in the manuscript are also found, including a model that
estimates a single correlation between the parallel coefficients in
the regressions of the two data types (rather than a single constant,
alpha; gi-plus-dbh-correlated.stan), and one that estimates a
correlation for each parallel coefficient
(gi-plus-dbh-multi-correlated.stan). We also considered analyzing
basal area increments rather than diameter/radial increments;
functions for reconstructing past basal area are found in the R driver
files, and models implementing analysis of basal area increments are
found in the branch basal-area-increments. Readers interested in these
alternative models are encouraged to contact M. E. K. Evans and/or
K. E. Holsinger.  Scripts for visualizing model output are found in
the repository as well, including autocorrelation.R, bivariate.R,
plot-mu-site.R, plot-results.R, and plot-site-effects.R. Additional
scripts for making figures are available from M. E. K. Evans.

