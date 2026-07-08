# bench_one.R <method> <N> [P]
# Fit ONE config on synthetic lattice ST data; print elapsed fit time (compile excluded
# where possible). Run as a subprocess under /usr/bin/time -l to capture peak RSS.
suppressMessages(library(rstan)); rstan_options(auto_write = TRUE)
set.seed(1)
args   <- commandArgs(trailingOnly = TRUE)
method <- args[1]                     # "dense" | "sparse"
N      <- as.integer(args[2])
P      <- ifelse(length(args) >= 3, as.integer(args[3]), 20L)
root   <- "/Users/licheng/Documents/DiSTect"
source(file.path(root, "repro/m0/neighbors.R"))

## --- synthetic lattice: square-ish grid, P genes, spatially-clustered disease ---
side <- ceiling(sqrt(N))
gx <- ((seq_len(N)-1) %%  side) + 1L
gy <- ((seq_len(N)-1) %/% side) + 1L
coords <- data.frame(x = gx, y = gy)
X <- matrix(rnorm(N*P), N, P)
beta_true <- c(1.5, -1.2, rep(0, P-2))            # 2 signal genes
# disease: smooth spatial field (via block pattern) + gene signal -> identifiable eta
field <- as.integer((gx %/% 5 + gy %/% 5) %% 2)   # checkerboard-ish spatial structure
mu <- X %*% beta_true + 1.5*field - 1
y  <- rbinom(N, 1, plogis(mu))
matrix_x <- cbind(X, x = coords$x, y = coords$y)

DENSE <- "
data { int<lower=0> N; int<lower=0> P; matrix[N,P+2] x; int y[N]; }
parameters { vector[P] beta; real<lower=0,upper=8> eta; vector<lower=0>[P] beta_gamma; real<lower=0,upper=1> w; }
model {
  vector[N] mu; matrix[N,N] prob_neigh;
  for(i in 1:P){ w~uniform(0,1); beta_gamma[i]~inv_gamma(5,50);
    target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]), log(w)+normal_lpdf(beta[i]|0,beta_gamma[i])); }
  for(i in 1:N){ for(j in 1:N){ eta~uniform(0,8);
    if (j!=i && sqrt(square(x[i,P+1]-x[j,P+1])+square(x[i,P+2]-x[j,P+2]))<=1) prob_neigh[i,j]=eta*y[j]; else prob_neigh[i,j]=0; }
    mu[i]=dot_product(x[i,1:P],beta)+sum(prob_neigh[i,]); }
  y ~ bernoulli_logit(mu);
}"
SPARSE <- "
data { int<lower=0> N; int<lower=0> P; matrix[N,P] x; int y[N]; vector[N] c; }
parameters { vector[P] beta; real<lower=0,upper=8> eta; vector<lower=0>[P] beta_gamma; real<lower=0,upper=1> w; }
model {
  vector[N] mu; w~uniform(0,1); eta~uniform(0,8);
  for(i in 1:P){ beta_gamma[i]~inv_gamma(5,50);
    target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]), log(w)+normal_lpdf(beta[i]|0,beta_gamma[i])); }
  mu = x*beta + eta*c; y ~ bernoulli_logit(mu);
}"

if (method == "dense") {
  m <- stan_model(model_code = DENSE)
  data <- list(N = N, P = P, x = matrix_x, y = y)
} else {
  m <- stan_model(model_code = SPARSE)
  cvec <- neighbor_sum(coords, y)
  data <- list(N = N, P = P, x = X, y = y, c = cvec)
}
gc(reset = TRUE)
t <- system.time(fit <- vb(m, data = data, algorithm = "fullrank",
                           seed = 128, iter = 3000, tol_rel_obj = 1e-5))["elapsed"]
cat(sprintf("RESULT method=%s N=%d P=%d fit_seconds=%.2f\n", method, N, P, t))
