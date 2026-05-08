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

#We divide by 100 to get the correct units in the example
env_array <-  example_2$env_array/100
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

    ##           BIO01     BIO10      BIO11     BIO12     BIO16      BIO17
    ## 2.5%  0.1221381 0.2139449 0.03847143 0.8808287 0.1014201 0.02213420
    ## 25%   0.1687778 0.2576851 0.08884324 1.1444185 0.1434298 0.04366665
    ## 50%   0.1871020 0.2663721 0.11528053 1.3157598 0.1740133 0.05776966
    ## 75%   0.2042655 0.2733272 0.14381416 1.5115700 0.2091240 0.07222351
    ## 97.5% 0.2429674 0.2829841 0.20728660 1.9581223 0.2955371 0.10721530

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
    ## [1,] 2421.222 2475.048 2483.906 2550.385 2580.205 2586.134 2593.560 2594.927
    ## [2,] 2368.020 2421.846 2430.704 2520.828 2550.648 2556.578 2540.358 2541.725
    ##          [,9]    [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 2598.267 2604.629 2605.736 2677.540 2688.451 2915.417 2945.273
    ## [2,] 2545.065 2551.427 2552.534 2647.984 2635.249 2885.860 2915.716

``` r

plot(model_BICs,model_AICs,type="p",xlab="BIC",ylab="AIC")
```

![](unuseable_models_files/figure-html/rank_by_AIC-1.png)

``` r

order(model_BICs)
```

    ##  [1] 11  8 14  5  3  1 13 15  7  9 10  2 12  4  6

``` r

order(model_AICs)
```

    ##  [1] 11  8 14  5 13 15  7  3  9 10  1 12  2  4  6

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

    ## [1] 1175.01

``` r

min(sapply(X=model_11_results, FUN=function(y){y$value}))
```

    ## [1] 1175.495

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

    ##      bestlogliks convergences parms_dists       mu1       mu2    sigltil1
    ## [1,]    1175.010            0     0.00000 0.2904673 0.3931544 0.002967830
    ## [2,]    1175.589            0    56.97587 0.2836055 0.3340812 0.003534804
    ## [3,]    1175.684            0    28.57079 0.2876660 0.3513049 0.002763906
    ## [4,]    1175.861            0    47.27432 0.2900079 0.3558192 0.002667892
    ## [5,]    1176.959            0    71.45096 0.2856940 0.3029788 0.002470630
    ## [6,]    1177.181            0   104.11978 0.2931011 0.3578675 0.002299245
    ## [7,]    1178.053            0   139.39522 0.2838179 0.2741103 0.002125990
    ## [8,]    1178.444            0   200.40759 0.2823220 0.2645135 0.001881315
    ##         sigltil2    sigrtil1   sigrtil2      ctil        pd     o_mat1
    ## [1,]  0.11063183 0.006786097 0.05675492 -9.639775 0.5532325 -0.9957053
    ## [2,]  0.37929728 0.006102938 0.04884209 -7.613571 0.5572592 -0.9953265
    ## [3,]  5.63506911 0.007317239 0.05301459 -7.698887 0.5491604 -0.9957239
    ## [4,]  0.05362719 0.007628213 0.02597972 -7.830341 0.5577326 -0.9950168
    ## [5,]  3.67097009 0.007848931 0.04500201 -6.317768 0.5634358 -0.9948601
    ## [6,]  0.12895618 0.008911315 0.05235808 -8.163370 0.5503899 -0.9945942
    ## [7,]  0.11473942 0.009266814 0.04025677 -5.053533 0.5685332 -0.9957438
    ## [8,] 38.77753040 0.009862485 0.03782166 -4.797914 0.5656953 -0.9970969
    ##          o_mat2      o_mat3     o_mat4
    ## [1,] 0.09257892 -0.09257892 -0.9957053
    ## [2,] 0.09656640 -0.09656640 -0.9953265
    ## [3,] 0.09237940 -0.09237940 -0.9957239
    ## [4,] 0.09970787  0.09970787  0.9950168
    ## [5,] 0.10125922 -0.10125922 -0.9948601
    ## [6,] 0.10383775 -0.10383775 -0.9945942
    ## [7,] 0.09216396 -0.09216396 -0.9957438
    ## [8,] 0.07614357 -0.07614357 -0.9970969

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

    ##      bestlogliks convergences parms_dists       mu1       mu2    sigltil1
    ## [1,]    1175.393            0     0.00000 0.2883993 0.3637573 0.003343087
    ## [2,]    1175.575            0    30.83822 0.2875756 0.3462229 0.003060580
    ## [3,]    1175.920            0    70.75508 0.2896679 0.3661489 0.002721052
    ## [4,]    1176.679            0   102.13801 0.2842533 0.2994654 0.002523407
    ## [5,]    1176.754            0   153.74904 0.2920729 0.3836868 0.002234328
    ## [6,]    1177.282            0    30.99508 0.3051778 0.4891399 0.003078746
    ## [7,]    1177.499            0   197.20093 0.2879086 0.3153094 0.002044076
    ## [8,]    1177.931            0   205.10175 0.2832989 0.2773803 0.002016125
    ##      sigltil2    sigrtil1   sigrtil2       ctil        pd     o_mat1     o_mat2
    ## [1,]      Inf 0.006260700 0.05263638  -8.852029 0.5465265 -0.9946534 0.10326924
    ## [2,]      Inf 0.006847563 0.05053014  -8.063806 0.5546072 -0.9948292 0.10156256
    ## [3,]      Inf 0.007064089 0.05150642  -9.377373 0.5356041 -0.9953305 0.09652530
    ## [4,]      Inf 0.007779143 0.04457290  -6.153619 0.5604980 -0.9954786 0.09498622
    ## [5,]      Inf 0.008355949 0.05346188  -9.851824 0.5413305 -0.9960327 0.08898753
    ## [6,]      Inf 0.006968360 0.07663110 -10.455178 0.5700551 -0.9939423 0.10990344
    ## [7,]      Inf 0.009314397 0.04832917  -6.112727 0.5648731 -0.9957713 0.09186720
    ## [8,]      Inf 0.009744197 0.04065154  -5.044924 0.5718964 -0.9967900 0.08006079
    ##           o_mat3     o_mat4
    ## [1,] -0.10326924 -0.9946534
    ## [2,] -0.10156256 -0.9948292
    ## [3,] -0.09652530 -0.9953305
    ## [4,] -0.09498622 -0.9954786
    ## [5,] -0.08898753 -0.9960327
    ## [6,] -0.10990344 -0.9939423
    ## [7,] -0.09186720 -0.9957713
    ## [8,] -0.08006079 -0.9967900

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

    ##         mu1       mu2  sigltil1  sigrtil1  sigrtil2       ctil        pd
    ## 1 0.2791391 0.2637573 -5.828900 -4.916322 -3.308011  -5.313797 0.2597900
    ## 2 0.2824500 0.3137573 -5.762449 -5.001852 -3.092502  -6.985897 0.2129451
    ## 3 0.2883993 0.3637573 -5.700861 -5.073463 -2.944348  -8.852029 0.1866459
    ## 4 0.2903715 0.4137573 -5.734683 -5.050950 -2.823565 -10.334806 0.1821756
    ## 5 0.2944830 0.4637573 -5.730928 -5.060225 -2.728970 -12.025495 0.1757218
    ## 6 0.2986233 0.5137573 -5.728509 -5.066494 -2.649442 -13.722028 0.1714030
    ##      o_par1 sigltil2
    ## 1 -3.235374      Inf
    ## 2 -3.233343      Inf
    ## 3 -3.245046      Inf
    ## 4 -3.230395      Inf
    ## 5 -3.229563      Inf
    ## 6 -3.228966      Inf

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

    ## [1] 0.2883993 0.3637573

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

    ##      bestlogliks convergences parms_dists        mu    sigltil    sigrtil
    ## [1,]    1255.414            0   0.0000000 0.2126061 0.02382594  8.6525472
    ## [2,]    1255.415            0   0.1346771 0.2123234 0.02377837 86.8102291
    ## [3,]    1255.415            0   0.3871227 0.2123189 0.02377651  2.0310168
    ## [4,]    1255.415            0   0.5669172 0.2123105 0.02377739  1.4799466
    ## [5,]    1255.416            0   0.8378647 0.2123271 0.02378026  1.0533136
    ## [6,]    1255.416            0   0.5536659 0.2122974 0.02380050  1.4997267
    ## [7,]    1255.416            0   0.9574630 0.2123270 0.02377991  0.9350686
    ## [8,]    1255.416            0   1.0780076 0.2123658 0.02378607  0.8395000
    ##           ctil        pd o_mat
    ## [1,] -1.785504 0.6404550     1
    ## [2,] -1.769444 0.6422606     1
    ## [3,] -1.768965 0.6423930     1
    ## [4,] -1.768070 0.6425047     1
    ## [5,] -1.769078 0.6424098     1
    ## [6,] -1.759311 0.6433603     1
    ## [7,] -1.769112 0.6424351     1
    ## [8,] -1.771164 0.6422635     1

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

    ##      bestlogliks convergences parms_dists        mu    sigltil sigrtil
    ## [1,]    1255.413            0 0.000000000 0.2131674 0.02396865     Inf
    ## [2,]    1255.413            0 0.009195052 0.2131345 0.02396340     Inf
    ## [3,]    1255.413            0 0.029997360 0.2131203 0.02395145     Inf
    ## [4,]    1255.413            0 0.043373134 0.2130857 0.02394376     Inf
    ## [5,]    1255.413            0 0.103083476 0.2129490 0.02390962     Inf
    ## [6,]    1255.413            0 0.026347380 0.2131802 0.02398354     Inf
    ## [7,]    1255.413            0 0.089208506 0.2129198 0.02391783     Inf
    ## [8,]    1255.413            0 0.118132855 0.2129092 0.02390103     Inf
    ##           ctil        pd o_mat
    ## [1,] -1.798828 0.6399018     1
    ## [2,] -1.797821 0.6397775     1
    ## [3,] -1.800396 0.6393561     1
    ## [4,] -1.799506 0.6393460     1
    ## [5,] -1.794903 0.6399583     1
    ## [6,] -1.794062 0.6408111     1
    ## [7,] -1.788930 0.6406384     1
    ## [8,] -1.793885 0.6398981     1

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
