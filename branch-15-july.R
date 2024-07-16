#### EM using logistic distribution ####

### Libraries ###
library(ggplot2)
library(quantreg)
library(splines)
library(conflicted)
library(rmetalog)
library(tidyr)
library(purrr)
library(dplyr)

### Function Declarations ###

# Function: l (likelihood) ------------------------------------------------

l <- function(par, w, X, y, scale) {
  # y_hat <- X %*% par
  
  # res   <- y - y_hat
  
  # s     <- sqrt(sum(res^2 * w) / sum(w))

  fs    <- dlogis(y, X %*% par, scale)
  
  wfs   <- w * log( fs )
  
  return(-sum(wfs))
}

debug_args <- function(...) {
  K       <<- 2 # components to estimate
  form    <<- y ~ t
  dat     <<- distr_1
  res_var <<- "y"
  id_col  <<- "id"
}

debug_args()

# Function start ----------------------------------------------------------

f_rhs <- formula(paste0("~", attr(terms(form),"term.labels"))) # get rhs for model matrix

## Initialize parameters

init_pars <-  get.init.pars(K = K, dat = dat, form = form, grp = "1", grp_vars = "t", id_col = "id")
# Get coefficients
betas <- init_pars$betas
scale <- init_pars$scale
X <- model.matrix(f_rhs, data = dat)
y <- dat[,res_var]

id_vector <- dat[,id_col]
tol <- 0.001
verbatim <- TRUE

## Init variables ##
f <- H <- P <- matrix(nrow = nrow(dat), ncol = K)

rep_measur <- ifelse(!is.null(id_vector), TRUE, FALSE) # repeated measurement indicator

if (rep_measur) {
  f_ind <- H <- matrix(nrow = length(unique(id_vector)), ncol = K)
}

Pi <- rep(1/K,K)
iter <- 0
fits <- vector("list",K)
ll <- -.Machine$double.xmax
conv <- 99

while (conv == 99) {
  
  
  ### Expectation ###
  
  iter <-  iter + 1

  ## Posterior Probabilities
  
  if (!rep_measur) {
    for (i in 1:K) {
      f[,i] <- dlogis(y, X %*% betas[,i], scale[i])
      H[,i] <- Pi[i] * f[,i]
    }
    P <- H / rowSums(H)
  } else {
    for (i in 1:K) {
      f[,i] <- dlogis(y, X %*% betas[,i], scale[i])
      f_ind[,i] <- tapply(f[,i], id_vector, prod) # likelihoods are multiplied
      H[,i] <- Pi[i] * f_ind[,i]
    }
    P <- H / rowSums(H)
    P <- P[id_vector,] # Expand P to original rank
  }

  ## Update Priors
  
  Pi <- colSums(P) / sum(P)

  ## Update scale parameter
 # res <- y - X %*% betas
  #scale <- sqrt(colSums(res^2 * P) / colSums(P))
  #names(scale) <- paste0("comp_",seq_len(K))

  ## Maximization
  for (i in seq_len(K)) {
    fits[[i]] <- optim(par = betas[,i], fn = l, X = X, y = y, w = P[,i], scale = scale[i])
  }
  ## Update parameters
  betas <- sapply(seq_len(K), \(i) fits[[i]]$par)
  colnames(betas) <- paste0("comp_",seq_len(K))
  
  old_ll <- ll
  ll <- sum(log(rowSums(H)))
  
  delta_ll <- ll - old_ll
  
  if (delta_ll < tol) conv <- 0
  
  if (verbatim) {
    cat("Iteration: ",iter," with Log-Likelihood: ", ll, "\n",sep = "")
  }
  
}

if (conv == 0) cat("Convergence after ", iter, " iterations. \n", sep = "")
if (conv == 1) cat("Failure to converge after ", iter, " iterations. \n", sep = "")


### Return object ###
model_obj <- list()

## Hard Classifications
comp_hat <- factor(paste0("comp_",apply(P, 1, \(x) which(x == max(x)))))
model_obj$comp_hat <- comp_hat 

model_obj$coefs <- betas

## Likelihood
model_obj$likelihood <- ll

## Original Data with estimated component
model_obj$data <- cbind(dat, comp_hat)

## Estimate conditional metalogs

df_metalogs <- model_obj$data |> 
  select(y, t, comp_hat) |> 
  nest(y = y) |> 
  mutate(ml_obj = map(y, \(y) metalog(y$y, term_limit = 5, term_lower_bound = 5)))

model_obj$cond_metalogs <- df_metalogs

# PREDICT -----------------------------------------------------------------

new_data <- data.frame(t = 1:4)

# Model matrix

# Extract rhs of formula for model matrix for predictions
X_new <- model.matrix(f_rhs, data = new_data)

# Predicted values
Y_hat <- X_new %*% model_obj$coefs

# Return data frame (in long format)
df_predicted <- cbind(new_data, Y_hat) |> 
  pivot_longer(matches("comp_"),
               names_to = "comp_hat",
               values_to = "y_hat")

df_densities <- model_obj$cond_metalogs |> 
  mutate(dens = map2(y, ml_obj, \(y, ml) dmetalog(ml, y$y, term = 5))) |> 
  select(-ml_obj) |> 
  unnest(c(y, dens))

# plot predicted distributions
ggplot() +
  geom_point(data = model_obj$data, mapping = aes(x = t, y = y, color = comp_hat)) +
  geom_ribbon(data = df_densities,
              mapping = aes(xmin = t, xmax = t + dens, y = y, group = interaction(comp_hat,t), fill = comp_hat), 
              alpha = 0.3)

