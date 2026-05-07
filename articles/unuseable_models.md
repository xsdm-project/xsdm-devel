# Unusable models: when a model does not have a maximum likelihood

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

#We divide by 100 to get the correct units in the example
env_array <-  example_2$env_array/100
dim(env_array)
```

    ## [1] 2728   39    6

``` r

dimnames(env_array)[[3]]
```

    ## NULL

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

    ##        [,1]  [,2]  [,3]   [,4]  [,5]  [,6]
    ## 2.5%  12.21 21.39  3.84  88.08 10.14  2.21
    ## 25%   16.87 25.76  8.88 114.44 14.34  4.36
    ## 50%   18.71 26.63 11.52 131.57 17.40  5.77
    ## 75%   20.42 27.33 14.38 151.15 20.91  7.22
    ## 97.5% 24.29 28.29 20.72 195.81 29.55 10.72

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
    ## [1,] 2419.798 2473.469 2482.653 2550.380 2580.196 2586.113 2593.530 2595.866
    ## [2,] 2366.596 2420.267 2429.451 2520.824 2550.639 2556.556 2540.328 2542.664
    ##          [,9]    [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 2596.864 2602.629 2604.562 2676.690 2683.136 2915.412 2945.931
    ## [2,] 2543.662 2549.427 2551.360 2647.133 2629.934 2885.855 2916.374

``` r

plot(model_BICs,model_AICs,type="p",xlab="BIC",ylab="AIC")
```

![](unuseable_models_files/figure-html/rank_by_AIC-1.png)

``` r

order(model_BICs)
```

    ##  [1] 11  8 14  5  3  1 13 15  7 10  9  2 12  4  6

``` r

order(model_AICs)
```

    ##  [1] 11  8 14  5 13 15  7 10  3  9  1 12  2  4  6

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

    ## [1] 1174.298

``` r

min(sapply(X=model_11_results, FUN=function(y){y$value}))
```

    ## [1] 1174.292

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

    ##      bestlogliks convergences parms_dists      mu1      mu2  sigltil1 sigltil2
    ## [1,]    1174.298            1    0.000000 29.73806 52.15672 0.6490663 7.125586
    ## [2,]    1174.483            1    7.225786 29.28605 45.31537 0.6433075 6.385971
    ## [3,]    1174.496            0    8.340648 29.21119 44.27018 0.6421165 6.270425
    ## [4,]    1174.544            0    9.543757 29.13078 43.12185 0.6373497 6.129187
    ## [5,]    1174.550            1    9.351745 29.13979 43.29941 0.6436414 6.143396
    ## [6,]    1174.556            1    9.486022 29.13161 43.17355 0.6400864 6.128443
    ## [7,]    1174.699            1   11.492589 28.96350 41.21234 0.6330165 5.833009
    ## [8,]    1174.742            0   13.247792 28.85798 39.63815 0.6447839 5.712863
    ##       sigrtil1    sigrtil2      ctil        pd    o_mat1      o_mat2     o_mat3
    ## [1,] 0.3141714    164.7586 -14.00756 0.5419413 0.9967329 -0.08076818 0.08076818
    ## [2,] 0.3212946   1102.1142 -11.72764 0.5432858 0.9963726 -0.08509729 0.08509729
    ## [3,] 0.3204658 304712.1682 -11.34549 0.5435118 0.9963343 -0.08554551 0.08554551
    ## [4,] 0.3228786    650.1384 -10.99477 0.5436531 0.9962384 -0.08665501 0.08665501
    ## [5,] 0.3210504    134.1798 -11.06821 0.5439611 0.9962915 -0.08604206 0.08604206
    ## [6,] 0.3251303  47722.9586 -11.02311 0.5434827 0.9962437 -0.08659421 0.08659421
    ## [7,] 0.3244256    194.2317 -10.58886 0.5408163 0.9962035 -0.08705517 0.08705517
    ## [8,] 0.3214258    929.1494  -9.76396 0.5458320 0.9961687 -0.08745205 0.08745205
    ##         o_mat4
    ## [1,] 0.9967329
    ## [2,] 0.9963726
    ## [3,] 0.9963343
    ## [4,] 0.9962384
    ## [5,] 0.9962915
    ## [6,] 0.9962437
    ## [7,] 0.9962035
    ## [8,] 0.9961687

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
    ## [1,]    1174.331            0   0.0000000 29.63142 49.81497 0.3241713      Inf
    ## [2,]    1174.334            0   0.5133777 29.58757 49.47035 0.3155704      Inf
    ## [3,]    1174.341            0   1.3445752 29.56323 48.68139 0.3167692      Inf
    ## [4,]    1174.348            0   0.1555118 29.60017 49.93797 0.3228371      Inf
    ## [5,]    1174.394            0   3.3201461 29.48731 46.81444 0.3145308      Inf
    ## [6,]    1174.405            0   2.5701106 29.39857 47.43588 0.3221710      Inf
    ## [7,]    1174.407            0   3.4990454 29.38082 46.60517 0.3238167      Inf
    ## [8,]    1174.422            0   2.4624018 29.42377 47.48744 0.3207519      Inf
    ##       sigrtil1 sigrtil2      ctil        pd     o_mat1     o_mat2      o_mat3
    ## [1,] 0.6261384 6.813339 -13.52742 0.5393367 -0.9964067 0.08469714 -0.08469714
    ## [2,] 0.6456272 6.833004 -13.16209 0.5416787 -0.9965494 0.08300215 -0.08300215
    ## [3,] 0.6455961 6.772182 -12.81278 0.5426067 -0.9964368 0.08434290 -0.08434290
    ## [4,] 0.6279639 6.811197 -13.61624 0.5386883 -0.9965230 0.08331848 -0.08331848
    ## [5,] 0.6500178 6.590210 -12.11781 0.5441018 -0.9962326 0.08672170 -0.08672170
    ## [6,] 0.6346573 6.582420 -12.58386 0.5409006 -0.9965286 0.08325072 -0.08325072
    ## [7,] 0.6357972 6.530596 -12.15736 0.5428998 -0.9963807 0.08500316 -0.08500316
    ## [8,] 0.6358889 6.547025 -12.75205 0.5404907 -0.9964804 0.08382619 -0.08382619
    ##          o_mat4
    ## [1,] -0.9964067
    ## [2,] -0.9965494
    ## [3,] -0.9964368
    ## [4,] -0.9965230
    ## [5,] -0.9962326
    ## [6,] -0.9965286
    ## [7,] -0.9963807
    ## [8,] -0.9964804

The best likelihoods obtained for this boundary model are very similar
to those obtained for the initial, non-boundary model, and the boundary
model has one fewer parameter, so we tentatively adopt the boundary
model over the earlier model. However, these optimization results reveal
that we probably still have not successfully optimized the likelihood,
or that the likelihood surface may have a ridge or an asymptote or other
pathological feature. One sees this by observing that whereas the top
several optimization results are similar in likelihood, they are spread
out in parameter space (the `parms\_dists` column shows distance in
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

![](unuseable_models_files/figure-html/plot_profiles-1.png)

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
    ## 1 29.52287 47.31497 -1.128312 -0.4544003 1.892310 -12.34334 0.1748812 28.18656
    ## 2 29.52699 47.36497 -1.128284 -0.4544711 1.893143 -12.36030 0.1748323 28.18656
    ## 3 29.53112 47.41497 -1.128253 -0.4545432 1.893976 -12.37725 0.1747863 28.18657
    ## 4 29.53525 47.46497 -1.128229 -0.4546104 1.894807 -12.39421 0.1747358 28.18657
    ## 5 29.53939 47.51497 -1.128196 -0.4546805 1.895637 -12.41115 0.1746875 28.18658
    ## 6 29.54351 47.56497 -1.128169 -0.4547517 1.896465 -12.42810 0.1746426 28.18659
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

![](unuseable_models_files/figure-html/make_pairs_plot_again-1.png)

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

![](unuseable_models_files/figure-html/look_at_inferred_ge_func-1.png)

``` r

param_list_bio$mu
```

    ## [1] 29.63142 49.81497

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
    ## [1,]    1255.412            0 0.0000000000 21.32529 2.399255 2.937853e+06
    ## [2,]    1255.412            0 0.0004470471 21.32563 2.399293 2.169292e+10
    ## [3,]    1255.412            0 0.0016974800 21.32659 2.399628 9.222548e+02
    ## [4,]    1255.412            0 0.0029701217 21.32678 2.399668 3.905636e+02
    ## [5,]    1255.412            0 0.0046965857 21.32545 2.399292 2.130360e+02
    ## [6,]    1255.412            0 0.0053307845 21.32748 2.399726 2.089303e+02
    ## [7,]    1255.412            0 0.0054598717 21.32548 2.399325 1.833251e+02
    ## [8,]    1255.412            0 0.0124665944 21.33611 2.401907 1.894166e+02
    ##           ctil        pd o_mat
    ## [1,] -1.806417 0.6390113     1
    ## [2,] -1.806707 0.6389671     1
    ## [3,] -1.806602 0.6390048     1
    ## [4,] -1.806654 0.6390159     1
    ## [5,] -1.806391 0.6390457     1
    ## [6,] -1.807263 0.6389515     1
    ## [7,] -1.806277 0.6390752     1
    ## [8,] -1.809611 0.6387623     1

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
    ## [1,]    1255.412            0 0.0000000000 21.32484 2.399158     Inf -1.806229
    ## [2,]    1255.412            0 0.0003132082 21.32509 2.399190     Inf -1.806409
    ## [3,]    1255.412            0 0.0001353714 21.32497 2.399183     Inf -1.806271
    ## [4,]    1255.412            0 0.0001059136 21.32477 2.399108     Inf -1.806307
    ## [5,]    1255.412            0 0.0005274700 21.32531 2.399242     Inf -1.806466
    ## [6,]    1255.412            0 0.0005119267 21.32530 2.399250     Inf -1.806441
    ## [7,]    1255.412            0 0.0005094554 21.32530 2.399253     Inf -1.806431
    ## [8,]    1255.412            0 0.0005185020 21.32531 2.399247     Inf -1.806447
    ##             pd o_mat
    ## [1,] 0.6390291     1
    ## [2,] 0.6390088     1
    ## [3,] 0.6390267     1
    ## [4,] 0.6390219     1
    ## [5,] 0.6390059     1
    ## [6,] 0.6390094     1
    ## [7,] 0.6390097     1
    ## [8,] 0.6390102     1

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

![](unuseable_models_files/figure-html/profile_model_5-1.png)

The profiles look adequate. So the model is suitable for use, unlike the
previous model. The BIC and AIC for this model are valid. If the second-
or third-best model is suitable in the same manner, the lowest-BIC model
suitable model should tend to be preferred.
