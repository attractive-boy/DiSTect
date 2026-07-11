# Draft email: DiSTect btaf530 simulation reproduction query

To: qihuang.zhang@mcgill.ca

Subject: Reproducing DiSTect (btaf530) Simulation 1 -- a question about the eta standard-error columns in Tables C1-C3

Dear Prof. Zhang,

I have been reproducing the simulation study in Zhao, Deng, and Zhang, "DiSTect:
a Bayesian spatial model for disease-associated gene discovery and prediction in
spatial transcriptomics" (Bioinformatics, btaf530), using the public resources
(GitHub StaGill/DiSTect, Zenodo 17127211, and the supplementary PDF).

The eta **point estimates** reproduce well -- roughly unbiased on average, and
the proposed model clearly beats the naive model on signal-gene bias, as the
paper shows. My question is specifically about the eta **uncertainty** columns
(avgSEE / avgSEM) in Tables C1-C3, which I have not been able to reproduce or,
as far as I can tell, to achieve with any estimator on the stated design.

## The observation

For Simulation 1 (30x30 lattice, 20 covariates X ~ N(0,1),
beta = (1,2,3,-4,-5,0_15), 2000-sweep Gibbs, neighbor-average model), I estimate
the Cramer-Rao floor for the eta standard error to be about **0.11**, whereas
Table C1 (eta = 0.4) reports avgSEE = avgSEM = 0.044 (ADVI) and 0.021 (NUTS) --
2.5x to 5x below that floor. (The attached figure, `fig_eta_se_floor.png`, shows
the eta SE plateauing at ~0.11 as the covariate scale is reduced, well above the
reported values.)

I checked this four independent ways, all on your exact neighbor-average
convention and priors (Table C7: b1=5, b2=50, c1=8, v0=1e-6):

- GLM conditional pseudolikelihood (n=200): eta empirical SD ~0.32.
- Your spike-slab model via ADVI (both meanfield and fullrank, n=200): eta
  avgSEE ~0.18-0.25.
- The same model via NUTS with adapt_delta=0.99 (divergence-free): eta empirical
  SD ~0.24, mean posterior SD ~0.26 (n=8 datasets).
- A full-likelihood Monte-Carlo MLE (Geyer-Thompson) on the proper autologistic
  joint: its SE(eta) equals the pseudolikelihood SE (ratio ~1.0) across your
  three eta settings, i.e. the pseudolikelihood already attains the Cramer-Rao
  bound -- there is no more efficient estimator available.

Crucially, the floor does not depend on the covariate scale: sweeping the
covariate SD toward zero only lowers the eta SE to ~0.11, never to ~0.04.

## A related internal check

For the Table C1 Proposed (ADVI) eta row, avgBias = 0.087 is about twice
avgSEE = avgSEM = 0.044, yet avgCR = 95%. A roughly-unbiased estimator with SD
0.044 would need its bias to be well below 0.044 to reach 95% coverage; a bias
of ~2 SEs would push coverage toward ~50%. Achieving 95% coverage instead
implies a true SD near 0.13, which matches the ~0.11 floor above rather than the
reported 0.044. So the avgSEM and avgCR columns look mutually inconsistent to me.

## The same pattern in Simulations 2-3 and in the HER2 analysis

This is not specific to Simulation 1. The reported eta empirical SE is similarly
5x-27x below what I obtain in Simulations 2 and 3. It also appears in the
real-data analysis: the paper reports eta = 2.936 with a 95% interval
[2.926, 2.946] (width 0.02), whereas my reproduction on the same public HER2ST
data gives eta = 2.665 with interval [1.589, 3.999] (width ~2.4) -- roughly 120x
wider. This is why I would like to understand the SE / credible-interval
computation before interpreting the eta uncertainties; the point estimates
themselves are close.

## Questions

1. How were avgSEE and avgSEM defined/computed for the tables? In particular,
   are they the SD of the 200 point estimates and the mean posterior SD, or a
   Monte-Carlo standard error of the mean (SD/sqrt(200)) or similar? (SD/sqrt(200)
   from my ~0.24 empirical SD would be ~0.017, close to the reported 0.021 --
   which is why I suspect a definitional difference rather than an error.)
2. Was the Simulation-1 design exactly a single 30x30 lattice (900 spots) per
   replicate, or were more spots / multiple lattices pooled?
3. Were X generated as standard normal with beta = (1,2,3,-4,-5,...), or was the
   covariate / linear-predictor scale standardized in some way? (With
   X ~ N(0,1) the linear predictor has SD ~7.4, which strongly separates the
   logistic fit; I wondered whether an effective rescaling was applied.)
4. Were the data generated with the normalized neighbor average in Equation (2),
   or the unnormalized neighbor sum used in the released package
   (prob_neigh[i,j] = eta * y[j])?
5. Were Tables C1-C6 produced by scripts outside the public GitHub/Zenodo
   archive? If so, would you be willing to share them?

Separately, I am also working to reproduce the HER2-positive breast cancer and
STARmap PLUS (Alzheimer's) analyses, where two data points are blocking me:
(a) the disease-status labels actually used for the HER2ST sections A-H (the
public HER2ST tree carries pathologist annotations only for
A1/B1/C1/D1/E1/F1/G2/H1), and (b) access to the STARmap PLUS data (Zenodo record
5842625 and Single Cell Portal SCP1375 both appear to require authorization).
Any pointer to these would be a great help.

I would be glad to send the full reproduction bundle (scripts, logs, figures,
and the four diagnostics above) on request. Thank you very much for your time --
I want to make sure the discrepancy is on my side before drawing any conclusion.

Best regards,

Cheng Li
School of Informatics, Hunan University of Chinese Medicine
