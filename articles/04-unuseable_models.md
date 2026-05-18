# Unusable models: when a model does not have a maximum likelihood

\\ \newcommand{\mean}\[1\]{\overline{#1}} \newcommand{\var}{\text{var}}
\newcommand{\cov}{\text{cov}} \newcommand{\cor}{\text{cor}}
\newcommand{\Rp}{\text{Re}} \newcommand{\E}{\text{E}}
\newcommand{\ltsgr}{\text{ltsgr}} \newcommand{\expit}{\text{expit}}
\newcommand{\logit}{\text{logit}} \\

**Abstract.** When an xsdm model does not have a maximum likelihood, it
should not be used. This document shows an example of how that can
occur, and how to diagnose it. Along the way, another type of boundary
model (in addition to the \\p_d=1\\ case described elsewhere) is
illustrated, where sigma parameters are set to infinity, corresponding
to insensitivity of annual net growth to environmental changes in a
certain direction in environment space.

The Eastern glass lizard, *Ophisaurus ventralis*, is a legless lizard
found in the southeastern United States. It is the longest and heaviest
species of its genus, growing up to 108cm in total length.

We start by loading the data:

``` r

library(xsdm)
env_array <-  example_2$env_array

dim(env_array)
```

    ## [1] 2728   39    6

``` r

dimnames(env_array)[[3]]
```

    ## [1] "BIO01" "BIO10" "BIO11" "BIO12" "BIO16" "BIO17"

``` r

occ <-  example_2$occ_vec
length(occ)
```

    ## [1] 2728

Here, there are 6 environmental variables recorded for 39 years in 2728
locations, with accompanying detections and pseudo-absences in the
variable `occ`. The first three environmental variables (BIO1, BIO10,
BIO11) are temperature variables, and the last three (BIO12, BIO16,
BIO17) are precipitation variables. BIO1 is mean annual temperature,
BIO10 is mean temperature of the warmest quarter, BIO11 is mean
temperature of the coldest quarter, BIO12 is annual precipitation, BIO16
is precipitation of the wettest quarter, and BIO17 is precipitation of
the driest quarter.

Now look at the distributions of values of environmental variables to
make sure they are not on very different scales, which would cause
problems for optimization:

``` r

apply(FUN=quantile, X=env_array, MARGIN=3,prob=c(.025,.25,.5,.75,.975))
```

    ##          BIO01    BIO10     BIO11     BIO12    BIO16     BIO17
    ## 2.5%  12.21381 21.39449  3.847143  88.08287 10.14201  2.213420
    ## 25%   16.87778 25.76851  8.884324 114.44185 14.34298  4.366665
    ## 50%   18.71020 26.63721 11.528053 131.57598 17.40133  5.776966
    ## 75%   20.42655 27.33272 14.381416 151.15700 20.91240  7.222351
    ## 97.5% 24.29674 28.29841 20.728660 195.81223 29.55371 10.721530

These distributions look basically OK.

Now fit 15 models, each from 25 starting conditions:

``` r

models <- matrix(c(1,0,0,0,0,0,
                  0,1,0,0,0,0,
                  0,0,1,0,0,0,
                  0,0,0,1,0,0,
                  0,0,0,0,1,0,
                  0,0,0,0,0,1,
                  1,0,0,1,0,0,
                  1,0,0,0,1,0,
                  1,0,0,0,0,1,
                  0,1,0,1,0,0,
                  0,1,0,0,1,0,
                  0,1,0,0,0,1,
                  0,0,1,1,0,0,
                  0,0,1,0,1,0,
                  0,0,1,0,0,1), nrow=15, byrow=TRUE)
all_model_results <- list()
for (i in 1:nrow(models))
{
  print(i)
  env_dat <-  env_array[ , , models[i,]==1, drop = FALSE]
  starts <-  xsdm::start_parms(env_dat,num_starts=25)
  all_optim_results <-  list()
  for (j in 1:nrow(starts))
  {
    all_optim_results[[j]] <-  optim(par = starts[j,],
                                   fn = xsdm::loglik_math,
                                   method = "BFGS",
                                   env_dat = env_dat,
                                   occ = occ,
                                   negative=TRUE,
                                   control = list(trace=0)
                                   )
  }
  all_model_results[[i]] <- all_optim_results
}
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15

Within each model, rank the optimization results:

``` r

for (i in 1:length(all_model_results))
{
  values <-  sapply(X = all_model_results[[i]], FUN=function(x){x$value})
  inds <-  order(values)
  all_model_results[[i]] <-  all_model_results[[i]][inds]
}
```

Rank the models by BIC, bearing in mind that we’ve been working with the
negative of the likelihood:

``` r

model_BICs <- sapply(X=all_model_results,
                      FUN=function(x){
                        best_loglik = x[[1]]$value
                        num_parms = length(x[[1]]$par)
                        n = length(occ)
                        BIC = 2*best_loglik + num_parms*log(n)
                        return(BIC)
                      }
                    )
```

Also by AIC, then display:

``` r

model_AICs <- sapply(X=all_model_results,
                      FUN=function(x){
                        best_loglik = x[[1]]$value
                        num_parms = length(x[[1]]$par)
                        AIC = 2*best_loglik + 2*num_parms
                        return(AIC)
                      }
                    )

inds <- order(model_BICs)
rbind(model_BICs[inds],model_AICs[inds])
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
    ## [1,] 2419.839 2473.841 2482.646 2550.383 2580.205 2586.126 2593.545 2594.406
    ## [2,] 2366.637 2420.639 2429.444 2520.826 2550.648 2556.569 2540.343 2541.204
    ##          [,9]    [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 2595.994 2602.754 2603.486 2676.651 2683.060 2915.417 2945.961
    ## [2,] 2542.792 2549.553 2550.284 2647.095 2629.858 2885.860 2916.404

``` r

plot(model_BICs,model_AICs,type="p",xlab="BIC",ylab="AIC")
```

![](04-unuseable_models_files/figure-html/rank_by_AIC-1.png)

``` r

order(model_BICs)
```

    ##  [1] 11  8 14  5  3  1 13  7 15 10  9  2 12  4  6

``` r

order(model_AICs)
```

    ##  [1] 11  8 14  5 13  7 15 10  9  3  1 12  2  4  6

The AIC and BIC results are pretty well aligned, and the four best
models are the same:

``` r

models[order(model_BICs)[1:4],]
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]    0    1    0    0    1    0
    ## [2,]    1    0    0    0    1    0
    ## [3,]    0    0    1    0    1    0
    ## [4,]    0    0    0    0    1    0

These four models all use the same predictor, the fifth one, and then
some use various possible second predictors.

We emphasize that these BIC and AIC values may or may not be meaningful,
since B/AIC can only be computed when the likelihood has been
effectively maximized. We will elaborate below, but for now we accept
these are pseudo-BIC and pseudo-AIC values, subject to later validation
or rejection.

Among the models we fitted, there is one clear winner in pseudo-BIC (the
model with the lowest pseudo-BIC). So let’s investigate it further.
Start by optimizing it a bit harder to see if we can do any better.

``` r

i <- 11
env_dat <- env_array[,,models[i,]==1,drop=FALSE]
starts <- xsdm::start_parms(env_dat, num_starts = 100)
model_11_results <- list()
for (j in 1:nrow(starts))
{
  model_11_results[[j]] <- optim(par=starts[j,],fn=xsdm::loglik_math,
                                method="BFGS",
                                env_dat = env_dat,
                                occ = occ,negative=TRUE,
                                control = list(trace=0))
}
all_model_results[[11]][[1]]$value
```

    ## [1] 1174.318

``` r

min(sapply(X=model_11_results, FUN=function(y){y$value}))
```

    ## [1] 1174.359

We did slightly better.

Now move forward by looking at the results for this model, starting by
writing a convenience function for examine optimization results:

``` r

examine_optim_results <- function(optim_results,mask=NULL)
{
  #put optimization results in order from best to worst
  bestlogliks <- sapply(X=optim_results,FUN=function(x){x$value})
  inds <- order(bestlogliks)
  bestlogliks <- bestlogliks[inds]
  optim_results <- optim_results[inds]

  #model convergence
  convergences <- sapply(X=optim_results,FUN=function(x){x$convergence})

  #compute distances to the best result in parameter space
  best_parms_math <- optim_results[[1]]$par
  parms_dists_to_best <- lapply(
    X=optim_results,
    FUN=function(x){
      xsdm::dist_between_params(
        x$par,
        best_parms_math,
        mask=mask,
        give_closest_rep=TRUE)
    }
  )
  parms_dists <- sapply(X=parms_dists_to_best, FUN=function(x){x$distance})
  bestparms <- sapply(X=parms_dists_to_best, FUN=function(x){unlist(x$representative)})

  #put it all together
  return(rbind(bestlogliks,convergences,parms_dists,bestparms))
}

h <- examine_optim_results(all_model_results[[11]])
t(h[ ,1:8])
```

    ##      bestlogliks convergences parms_dists      mu1      mu2 sigltil1  sigltil2
    ## [1,]    1174.318            1    0.000000 29.94396 53.85131 7.254098 0.3229644
    ## [2,]    1174.403            0    6.634318 29.47977 47.65828 6.661796 0.3222862
    ## [3,]    1174.440            0    7.759930 29.39858 46.58588 6.535730 0.3207017
    ## [4,]    1174.443            1    4.791790 29.59701 49.29914 6.753602 0.3139233
    ## [5,]    1174.447            1    7.538554 29.39473 46.80321 6.566780 0.3215116
    ## [6,]    1174.483            1    9.261663 29.33182 45.19059 6.404816 0.3224701
    ## [7,]    1174.485            1    8.926674 29.28884 45.47378 6.400157 0.3251839
    ## [8,]    1174.493            1    8.936779 29.30083 45.47238 6.411428 0.3211931
    ##          sigrtil1  sigrtil2      ctil        pd     o_mat1    o_mat2     o_mat3
    ## [1,] 1.175588e+12 0.6317486 -14.79527 0.5402512 0.08385240 0.9964782 -0.9964782
    ## [2,] 1.702684e+02 0.6383495 -12.46191 0.5433223 0.08495379 0.9963849 -0.9963849
    ## [3,] 2.929599e+01 0.6401084 -12.12501 0.5434468 0.08511438 0.9963712 -0.9963712
    ## [4,] 1.035099e+02 0.6491168 -13.34314 0.5399802 0.08361673 0.9964980 -0.9964980
    ## [5,] 6.623826e+01 0.6407272 -12.17777 0.5430792 0.08438747 0.9964330 -0.9964330
    ## [6,] 1.885932e+02 0.6401751 -11.57117 0.5447481 0.08705572 0.9962034 -0.9962034
    ## [7,] 8.379447e+01 0.6338410 -11.78324 0.5428916 0.08517863 0.9963657 -0.9963657
    ## [8,] 9.438802e+01 0.6404986 -11.75474 0.5430017 0.08501808 0.9963794 -0.9963794
    ##          o_mat4
    ## [1,] 0.08385240
    ## [2,] 0.08495379
    ## [3,] 0.08511438
    ## [4,] 0.08361673
    ## [5,] 0.08438747
    ## [6,] 0.08705572
    ## [7,] 0.08517863
    ## [8,] 0.08501808

The very large values of `sigltil2` suggest the boundary model where
this parameter is set to `Inf`, corresponding to a direction in
environment space along which annual net growth is insensitive to
changes in the environment.

So we consider the corresponding boundary model:

``` r

i <- 11
env_dat <- env_array[ , , models[i,] == 1, drop=FALSE]
mask <- c(sigltil2 = Inf)
new_starts <- xsdm::start_parms(env_dat[occ == 1, , , drop=FALSE],
                               mask = mask,
                               num_starts = 100)

bdry_optim_results <- list()
for (j in 1:nrow(new_starts))
{
  bdry_optim_results[[j]] <- optim(par = new_starts[j,],
                                  fn = xsdm::loglik_math,
                                  method = "BFGS",
                                  env_dat = env_dat,
                                  occ = occ,
                                  mask = mask,
                                  negative = TRUE,
                                  control = list(trace=0, maxit=500))
}
```

Now look at these results:

``` r

h <- examine_optim_results(bdry_optim_results, mask = mask)
t(h[ ,1:8])
```

    ##      bestlogliks convergences parms_dists      mu1      mu2  sigltil1 sigltil2
    ## [1,]    1174.443            0   0.0000000 29.39769 46.71887 0.3213298      Inf
    ## [2,]    1174.471            0   0.7637072 29.35938 45.95669 0.3188537      Inf
    ## [3,]    1174.510            0   1.7685957 29.25876 44.99372 0.3203260      Inf
    ## [4,]    1174.587            0   4.0255684 29.15446 42.89161 0.3182872      Inf
    ## [5,]    1174.669            0   5.5140390 28.99943 41.45830 0.3234842      Inf
    ## [6,]    1174.812            0   8.0057813 28.82822 39.13701 0.3221507      Inf
    ## [7,]    1198.768            0  26.8749638 24.98625 22.77155 0.1193033      Inf
    ## [8,]    1198.769            0  26.8701888 24.98671 22.77528 0.1193456      Inf
    ##         sigrtil1 sigrtil2       ctil        pd     o_mat1      o_mat2
    ## [1,]   0.6424894 6.572514 -12.077186 0.5439307 -0.9964114  0.08464179
    ## [2,]   0.6385956 6.427518 -12.064210 0.5413249 -0.9963339  0.08554991
    ## [3,]   0.6383937 6.320893 -11.713544 0.5417975 -0.9963769  0.08504799
    ## [4,]   0.6469665 6.117436 -10.853675 0.5447301 -0.9961787  0.08733883
    ## [5,]   0.6370810 5.914408 -10.473585 0.5437911 -0.9962009  0.08708522
    ## [6,]   0.6458965 5.653954  -9.570628 0.5463650 -0.9961308  0.08788350
    ## [7,] 236.4614500 2.796797  -2.325270 0.6155689  0.9992977 -0.03747134
    ## [8,] 589.1254234 2.797183  -2.327854 0.6153615  0.9992968 -0.03749587
    ##           o_mat3     o_mat4
    ## [1,] -0.08464179 -0.9964114
    ## [2,] -0.08554991 -0.9963339
    ## [3,] -0.08504799 -0.9963769
    ## [4,] -0.08733883 -0.9961787
    ## [5,] -0.08708522 -0.9962009
    ## [6,] -0.08788350 -0.9961308
    ## [7,] -0.03747134 -0.9992977
    ## [8,] -0.03749587 -0.9992968

The best likelihoods obtained for this boundary model are very similar
to those obtained for the initial, non-boundary model, and the boundary
model has one fewer parameter, so we tentatively adopt the boundary
model over the earlier model. However, these optimization results reveal
that we probably still have not successfully optimized the likelihood,
or that the likelihood surface may have a ridge or an asymptote or other
pathological feature. One sees this by observing that whereas the top
several optimization results are similar in likelihood, they are spread
out in parameter space (the `parms_dists` column shows distance in
parameter space to the top-likelihood result).

To investigate the model further, we profile:

``` r

pnames <-  names(xsdm::make_mask_names(2))
pnames <-  pnames[!(pnames %in% names(mask))]

values <-  sapply(X=bdry_optim_results, FUN=function(x){x$value})
inds <-  order(values)
bdry_optim_results <-  bdry_optim_results[inds]

all_profiles <- list()
linc <-  rep(0.05, length(pnames))
rinc <-  rep(0.05, length(pnames))
linc[8] <-  0.001
rinc[8] <-  0.001
for (counter in 1:length(pnames))
{
  all_profiles[[counter]] <-  xsdm::profile_likelihood(
                              profile_parameter = pnames[counter],
                              increment_left =linc[counter],
                              increment_right = rinc[counter],
                              num_steps_left = 50,
                              num_steps_right = 50,
                              alpha = 0.95,
                              optim_param_vector = bdry_optim_results[[1]]$par,
                              env_dat=env_dat,
                              occ = occ,
                              mask = mask,
                              num_threads = 6
                            )
}
names(all_profiles) <-  pnames
```

Now plot these profiles:

``` r

plot_tool <- function(ap,index)
  {
  x <- ap[[index]]$profile$value_math
  y <- ap[[index]]$profile$loglik
  xlab <- names(ap)[index]
  thresh <- ap[[index]]$threshold
  plot(x,y,
       type="o",xlab=xlab,
       ylab="Log likelihood")
  lines(range(x),rep(thresh,2),type="l",
        lty="dashed",col="red")
}

par(mfrow=c(3,3))
plot_tool(all_profiles, 1)
plot_tool(all_profiles, 2)
plot_tool(all_profiles, 3)
plot_tool(all_profiles, 4)
plot_tool(all_profiles, 5)
plot_tool(all_profiles, 6)
plot_tool(all_profiles, 7)
plot_tool(all_profiles, 8)
```

![](04-unuseable_models_files/figure-html/plot_profiles-1.png)

These profiles are not dome-shaped, and have other idiosyncratic
features, confirming that we had not effectively optimized the
likelihood. The `mu2`, `mu2`, and `ctil` profile, in particular, show
problems.

We do a pairs plot based on the `mu2` profile to investigate further:

``` r

p <- all_profiles[[2]]$parameters
head(p)
```

    ##        mu1      mu2  sigltil1   sigrtil1 sigrtil2      ctil        pd   o_par1
    ## 1 29.27045 44.21887 -1.127135 -0.4515529 1.837563 -11.29533 0.1781613 3.053299
    ## 2 29.27457 44.26887 -1.127101 -0.4516401 1.838493 -11.31225 0.1780996 3.053307
    ## 3 29.27869 44.31887 -1.127068 -0.4517272 1.839422 -11.32916 0.1780382 3.053315
    ## 4 29.28281 44.36887 -1.127034 -0.4518143 1.840349 -11.34607 0.1779772 3.053323
    ## 5 29.28692 44.41887 -1.127000 -0.4519010 1.841274 -11.36298 0.1779162 3.053331
    ## 6 29.29104 44.46887 -1.126966 -0.4519871 1.842197 -11.37990 0.1778556 3.053339
    ##   sigltil2
    ## 1      Inf
    ## 2      Inf
    ## 3      Inf
    ## 4      Inf
    ## 5      Inf
    ## 6      Inf

``` r

p <- p[,1:8]
pairs(p)
```

![](04-unuseable_models_files/figure-html/make_pairs_plot_again-1.png)

This pairs plot helps us see the tradeoff going on between parameters.
As `mu2` is increased, other parameters, especially `ctil`, change
monotonically. These results suggest a ridge in the likelihood surface
that appears to rise asymptotically along some path in parameter space
for which `mu2` is increasing. We have not identified a maximum of the
likelihood function, and it looks as though there may not be one if the
increase is indeed asymptotic.

We look at what the growth-environment function looks like for the best
parameters we have found so far, maybe that will give some insight into
what is going wrong:

``` r

param_list <- bdry_optim_results[[1]]$par
param_list["sigltil2"] <- Inf
param_list <- param_list[names(xsdm::make_mask_names(2))]
param_list_bio <- xsdm::math_to_bio(param_list)
xsdm::interpret_parameters(param_list = param_list_bio,
                          plot_indices = c(1,2), env_dat=env_dat, occ = occ)
```

![](04-unuseable_models_files/figure-html/look_at_inferred_ge_func-1.png)

``` r

param_list_bio$mu
```

    ## [1] 29.39769 46.71887

Estimated values of `mu2` are outside the range of the environmental
data. The actual distribution of the species is across Florida and along
the coasts of the Southeastern United States, and the southern extent of
the species is very likely constrained, on the southern end of the
range, by the Atlantic Ocean and the Gulf of Mexico rather than by
temperature and precipitation. Ultimately the modeling probably fails
here for these reasons. A Bayesian approach to the same problem could
likely address some of the limitations, here, by setting appropriate
priors. In the frequentist setting, one must discard the model. AIC and
BIC values are invalid when the likelihood has not been adequately
optimized and when the likelihood surface in the vicinity of the optimum
is not approximately a dome.

Immediate next steps should include examination of the second- and
third-best models according to pseudo-BIC, to see if they have similar
problems. For brevity and because we are just trying to illustrate
statistical principles and workflows, we skip straight to examining the
fourth-best model, which may be simpler because it only uses one
environmental variable.

``` r

i <- 5
model_5_optim_results <-  all_model_results[[i]]
h <- examine_optim_results(model_5_optim_results)
t(h[,1:8])
```

    ##      bestlogliks convergences  parms_dists       mu  sigltil      sigrtil
    ## [1,]    1255.413            0 0.0000000000 21.33103 2.399389 2.866696e+05
    ## [2,]    1255.413            0 0.0002897991 21.33103 2.399479 7.663422e+07
    ## [3,]    1255.413            0 0.0006045789 21.33145 2.399515 2.315692e+03
    ## [4,]    1255.413            0 0.0027395298 21.32836 2.398674 1.141085e+05
    ## [5,]    1255.413            0 0.0006674455 21.33139 2.399592 2.265492e+03
    ## [6,]    1255.413            0 0.0012886117 21.33210 2.399704 1.432889e+03
    ## [7,]    1255.413            0 0.0157284885 21.34585 2.402722 4.342884e+02
    ## [8,]    1255.413            0 0.0039797575 21.33150 2.399587 2.531599e+02
    ##           ctil        pd o_mat
    ## [1,] -1.806775 0.6389829     1
    ## [2,] -1.806487 0.6390119     1
    ## [3,] -1.806843 0.6389766     1
    ## [4,] -1.806160 0.6390095     1
    ## [5,] -1.806428 0.6390204     1
    ## [6,] -1.806898 0.6389765     1
    ## [7,] -1.811445 0.6387384     1
    ## [8,] -1.806585 0.6390307     1

One can see the boundary model with `sigrtil` set to `Inf` should be
considered:

``` r

env_dat <- env_array[,,models[i,]==1,drop=FALSE]
mask <- c(sigrtil1 = Inf)
new_starts <- xsdm::start_parms(env_dat[occ==1,,,drop=FALSE],mask=mask,
                               num_starts=100)

bdry_optim_results5 <- list()
for (j in 1:nrow(new_starts))
{
  bdry_optim_results5[[j]] <- optim(par=new_starts[j,],fn=xsdm::loglik_math,
                                method="BFGS",
                                env_dat=env_dat,occ=occ,mask=mask,negative=TRUE,
                                control=list(trace=0,maxit=500))
}
```

Examine these results:

``` r

h <- examine_optim_results(bdry_optim_results5,mask=mask)
t(h[,1:8])
```

    ##      bestlogliks convergences  parms_dists       mu  sigltil sigrtil      ctil
    ## [1,]    1255.413            0 0.000000e+00 21.33070 2.399347     Inf -1.806559
    ## [2,]    1255.413            0 5.017857e-05 21.33070 2.399333     Inf -1.806609
    ## [3,]    1255.413            0 4.299752e-05 21.33074 2.399363     Inf -1.806553
    ## [4,]    1255.413            0 1.526335e-04 21.33084 2.399381     Inf -1.806608
    ## [5,]    1255.413            0 1.445972e-04 21.33084 2.399383     Inf -1.806597
    ## [6,]    1255.413            0 1.434946e-04 21.33084 2.399385     Inf -1.806589
    ## [7,]    1255.413            0 1.568588e-04 21.33085 2.399386     Inf -1.806601
    ## [8,]    1255.413            0 1.517942e-04 21.33085 2.399387     Inf -1.806595
    ##             pd o_mat
    ## [1,] 0.6390027     1
    ## [2,] 0.6389988     1
    ## [3,] 0.6390028     1
    ## [4,] 0.6389998     1
    ## [5,] 0.6390002     1
    ## [6,] 0.6390012     1
    ## [7,] 0.6390000     1
    ## [8,] 0.6390000     1

``` r

values <- sapply(X=bdry_optim_results5, FUN=function(x){x$value})
inds <- order(values)
bdry_optim_results5 <- bdry_optim_results5[inds]
```

This model looks like it probably optimized well.

Next do profiles:

``` r

all_profiles <- list()
pnames <- names(xsdm::make_mask_names(1))
pnames <- pnames[pnames!="sigrtil1"]
linc <- c(0.05,0.025,0.05,0.05)
rinc <- c(0.15,0.025,0.05,0.05)
for (counter in 1:length(pnames))
{
  all_profiles[[counter]] <-  xsdm::profile_likelihood(
                              profile_parameter = pnames[counter],
                              increment_left  = linc[counter],
                              increment_right = rinc[counter],
                              num_steps_left  = 50,
                              num_steps_right = 50,
                              alpha = 0.95,
                              optim_param_vector = bdry_optim_results5[[1]]$par,
                              env_dat = env_dat,
                              occ = occ,
                              mask = mask,
                              num_threads = 6
                            )
}
names(all_profiles) <- pnames

par(mfrow=c(2,3))
plot_tool(all_profiles,1)
plot_tool(all_profiles,2)
plot_tool(all_profiles,3)
plot_tool(all_profiles,4)
```

![](04-unuseable_models_files/figure-html/profile_model_5-1.png)

The profiles look adequate. So the model is suitable for use, unlike the
previous model. The BIC and AIC for this model are valid. If the second-
or third-best model is suitable in the same manner, the lowest-BIC model
suitable model should tend to be preferred.
