library(castor)
library(TESS)
library(slouch)
library(progress)
library(ggplot2)
library(dplyr)

## SETTINGS

set.seed(123)
ntrees <- 1 ## Just use one tree (can't change this atm)
nsim <- 50 ## Number of trait simulations
ntaxas <- c(25, 100, 250, 500) ## number of taxa in the tree
mvs <- c(0.0, 0.2, 0.4, 0.8, 1.6) ## measurement variance, i.e., trait += Normal(mean = 0, var = mv)
upper <- log(2)/0.1 ## upper bounds for the alpha-parameter. I.e., we only optimize for an alpha inbetween what corresponds to the range of a half-life of between 0.1 and infinity

## Simulate trees
trees <- list()
for (i in seq_along(ntaxas)){
  trees[[i]] <- tess.sim.taxa.age(ntrees, ntaxas[[i]], 1.0, 1.0, 0.5)[[1]]
}

## Simulate trait data
data <- list()
for (i in seq_along(trees)){
  l <- castor::simulate_bm_model(trees[[i]], sigma = 1, include_tips = TRUE, Nsimulations = nsim)$tip_states
  data[[i]] <- l
}


## Progress bar
pb <- progress_bar$new(total = nsim * length(ntaxas) * length(mvs))

l1 <- list()
for (k in seq_along(mvs)){
  mv <- mvs[k]
  l2 <- list()
  for (i in seq_along(ntaxas)){
    ntaxa <- ntaxas[[i]]
    tree <- trees[[i]]
    l3 <- list()
    for (j in 1:nsim){
      d1 <- data[[i]][j,] + rnorm(ntaxa, mean = 0, sd = sqrt(mv))
      dat <- data.frame("trait" = d1)
      rownames(dat) <- tree$tip.label

      ## use phylolm because it is much more efficiently implemented than slouch.
      ## phylolm does not treat the mv as a fixed input, but instead models it as an unknown parameter in the model
      ## anyway, it shows some of the effects that we want to illustrate
      m0 <- phylolm(trait ~ 1,
                    phy = tree,
                    data = dat,
                    model = "OUfixedRoot",
                    measurement_error = TRUE,
                    upper.bound = upper)
      out0 <- c("ntaxa" = paste0("ntaxa = ", ntaxa), "mv" = paste0("mv = ", mv), "a" = m0$optpar, "mv_accounted" = "yes")

      m1 <- phylolm(trait ~ 1,
                    phy = tree,
                    data = dat,
                    model = "OUfixedRoot",
                    upper.bound = upper)
      out1 <- c("ntaxa" = paste0("ntaxa = ", ntaxa), "mv" = paste0("mv = ", mv), "a" = m1$optpar, "mv_accounted" = "no", upper.bound = upper)

      l3[[j]] <- bind_rows(out0, out1)
      pb$tick() ## update progressbar
    }
    l2[[i]] <- bind_rows(l3)
  }
  l1[[k]] <- bind_rows(l2)
}

## Clean up the dataframe
df <- bind_rows(l1)
df$a <- as.numeric(df$a)

df$ntaxa <- factor(df$ntaxa)
levels(df$ntaxa) <- paste0("ntaxa = ", ntaxas)

## Do the plot
yformat <- function(x) sprintf("%.1f", abs(round(log(2)/x, 2)))
p <- ggplot(df, aes(x = mv_accounted, fill = mv_accounted, y = a)) +
  geom_violin() +
  scale_y_sqrt(breaks = c(0.01980421, 0.2310491, 0.6931472, 1.732868, log(2)/0.2, 7), labels = yformat) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_classic() +
  facet_grid(cols = vars(ntaxa), rows = vars(mv)) +
  labs(x = "Observation variance accounted for in model", y = "Phylogenetic half-life (ML estimate)") +
  #coord_cartesian(ylim = c(0.0, 8.0)) +
  scale_fill_manual(values = c("black", "gray")) +
  theme(legend.position = "none")
p
ggsave(filename = "figures/ou_inference_bm_sim_measurement_errors.pdf", p, width = 200, height = 180, units = "mm")

# Draft of the paragraph

# As we mention in the main text, the effects of measurement error in comparative studies has been extensively studied.
# Several software packages have implemented functionality to incorporate measurement error as part of the model, and several of these (Cite ...) were available at the time when the study by Cooper et al. 2016 was presented.
# We will provide here a small example where we simulate trait data under various conditions, and try to recover the parameter estimates with and without accounting for measurement error in the model.
# Even though the effects of measurement error is known in the theoretical literature, we hope that this simulation example will be pedagogical for the readers who are not as acquainted with the more theoretical statistical literature.

# We begin by simulating a set of trees under the reconstructed birth-death model, with a speciation rate of 1.0, and an extinction rate of 0.5.
# We condition the tree simulations to result in a tree whose height is 1.0 time units, and for the taxa to be some fixed number, using the R-package TESS (Höhna et al 2015?).
# We simulate four trees of different size, by varying the conditioning to result in a tree with 25, 100, 250 and 500 taxa (with no replicates per setting).
#
# For each of the different trees, we use the R-package `castor` (cite Louca & Doebeli 2017) to simulate trait data under a Brownian motion with unit variance, with 50 replicates.
# From the simulated trait data we generate datasets that include measurement error with a variance that ranges from zero and upwards.
# In other words, we add a secondary source of variation on top of the simulated data.
# I.e., we let $y_i = y_{i,BM} + Normal(0, \sigma^2)$ for each species $i$ ($\sigma^2$ ranging from 0.0 to 1.6).
# Now, the trait value $y_i$ is composed of two elements: the phylogenetically structured BM data, and a noise component which is independently distributed with respect to the phylogeny.
# This is analogous to a scenario with real data, where we model the trait using some phylogenetically informed model, but we assume that there is some non-zero observation uncertainty.

# For each instance of the simulated BM data, we attempt to recover the parameter estimate of alpha by fitting a single-optimum Ornstein-Uhlenbeck model.
# To do so, we use the R-package `phylolm` (Ho & Ané 2014). Since the true value is $\alpha = 0$, we fit the model and plot the results in terms of $\alpha$.
# We set the limits of the numerical optimizer between $0 < $\alpha$ < \log(2)/0.1$, which corresponds to the range of half-life between 10% of the tree height and infinity.

## CAPTION

# Maximum-likelihood estimates of the phylogenetic half-life ($t_{1/2} = \log(2)/\alpha$) based on simulated data (50 replicates) under a range of conditions. A half-life of around 3 or greater represents strong phylogenetic signal (BM). A half-life of 0.1, of 10% of the tree height represents weak phylogenetic signal (OU). The y-axis is scaled such that $\sqrt{\alpha}$ is spaced linearly. The dashed horizontal line represents $t_{1/2} = \infty$, or $\alpha = 0$. The simulation and inference conditions include varying tree size (25, 100, 250, 500), varying measurement error added to the simulated trait data (Gaussian noise with variance equal to 0.0, 0.2, 0.4, 0.8 and 1.6), as well as whether the phylogenetic comparative method accounted for measurement error or not. We used the package R-package `phylolm`(Ho & Ané 2014) to fit the models. In order to account for observational error, `phylolm` adds a random variable $\sigma_\text{error}^2$ to the model. It will fit the model using generalized least-squares, a standard estimation technique in comparative phylogenetic methods, however it will partition the among-species covariance as $V_\text{total} = V_\text{OU} + \sigma_\text{error}^2 I$, where $V_\text{OU}$ is the covariance structure of the OU model, and $\sigma_\text{error}^2 I$ is a diagonal matrix representing the observation variance.



