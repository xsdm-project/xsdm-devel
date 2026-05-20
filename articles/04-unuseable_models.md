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
illustrated. In the new boundary model, sigma parameters are set to
infinity, corresponding to insensitivity of annual net growth to
environmental changes in a certain direction in environment space. This
example is based on occurrence data from GBIF for *Ophisaurus
ventralis*, the Eastern glass lizard.

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
make sure they are not on very different scales, which could cause
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
  env_dat <-  env_array[ , , models[i,]==1, drop = FALSE]
  starts <-  start_parms(env_dat[occ==1,,,drop=FALSE],num_starts=25)
  all_optim_results <-  list()
  for (j in 1:nrow(starts))
  {
    all_optim_results[[j]] <-  optim(par = starts[j,],
                                   fn = loglik_math,
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
    ## [1,] 2419.989 2473.809 2482.714 2550.383 2580.204 2586.123 2592.968 2593.547
    ## [2,] 2366.787 2420.607 2429.512 2520.826 2550.648 2556.567 2539.766 2540.345
    ##          [,9]    [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 2595.855 2602.753 2603.377 2676.652 2682.989 2915.417 2945.925
    ## [2,] 2542.653 2549.551 2550.175 2647.096 2629.787 2885.860 2916.368

``` r

plot(model_BICs,model_AICs,type="p",xlab="BIC",ylab="AIC")
```

![](04-unuseable_models_files/figure-html/rank_by_AIC-1.png)

``` r

order(model_BICs)
```

    ##  [1] 11  8 14  5  3  1  7 13 15 10  9  2 12  4  6

``` r

order(model_AICs)
```

    ##  [1] 11  8 14  5  7 13 15 10  9  3  1 12  2  4  6

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
starts <- start_parms(env_dat[occ==1,,,drop=FALSE], num_starts = 100)
model_11_results <- list()
for (j in 1:nrow(starts))
{
  model_11_results[[j]] <- optim(par=starts[j,],fn=loglik_math,
                                method="BFGS",
                                env_dat = env_dat,
                                occ = occ,negative=TRUE,
                                control = list(trace=0))
}
all_model_results[[11]][[1]]$value
```

    ## [1] 1174.393

``` r

min(sapply(X=model_11_results, FUN=function(y){y$value}))
```

    ## [1] 1174.361

About the same.

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
      dist_between_params(
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
    ## [1,]    1174.393            0   0.0000000 29.51105 48.25513 0.6423693 6.723245
    ## [2,]    1174.398            1   0.5760184 29.51613 48.71516 0.6354772 6.722903
    ## [3,]    1174.422            0   0.8199987 29.43727 47.45322 0.6372517 6.603519
    ## [4,]    1174.429            0   0.5052323 29.48902 47.79219 0.6233952 6.590016
    ## [5,]    1174.533            0   4.3878324 29.22367 44.09852 0.6418905 6.249746
    ## [6,]    1174.974            1  10.0403875 28.80064 38.67883 0.6453429 5.506687
    ## [7,]    1175.021            1  12.0833898 28.61807 36.80735 0.6437729 5.305104
    ## [8,]    1175.156            1  13.5499744 28.53063 35.44734 0.6559985 5.142161
    ##       sigrtil1   sigrtil2       ctil        pd    o_mat1      o_mat2     o_mat3
    ## [1,] 0.3195678  104.59916 -12.674704 0.5425995 0.9964647 -0.08401277 0.08401277
    ## [2,] 0.3224764 4135.29690 -13.019628 0.5408669 0.9965094 -0.08348103 0.08348103
    ## [3,] 0.3210268  390.88295 -12.521493 0.5419016 0.9964559 -0.08411682 0.08411682
    ## [4,] 0.3260189  238.90338 -12.859989 0.5396807 0.9963178 -0.08573728 0.08573728
    ## [5,] 0.3203065   35.25965 -11.299033 0.5434350 0.9962621 -0.08638192 0.08638192
    ## [6,] 0.3188111  434.20979  -9.742483 0.5443215 0.9960838 -0.08841417 0.08841417
    ## [7,] 0.3219691  224.83701  -8.912121 0.5466605 0.9961405 -0.08777286 0.08777286
    ## [8,] 0.3176739  599.86716  -8.362218 0.5483841 0.9961031 -0.08819685 0.08819685
    ##         o_mat4
    ## [1,] 0.9964647
    ## [2,] 0.9965094
    ## [3,] 0.9964559
    ## [4,] 0.9963178
    ## [5,] 0.9962621
    ## [6,] 0.9960838
    ## [7,] 0.9961405
    ## [8,] 0.9961031

The very large values of `sigrtil1` suggest the boundary model where
this parameter is set to `Inf`, corresponding to a direction in
environment space along which annual net growth is insensitive to
changes in the environment.

So we consider the corresponding boundary model:

``` r

i <- 11
env_dat <- env_array[ , , models[i,] == 1, drop=FALSE]
mask <- c(sigrtil1 = Inf)
new_starts <- start_parms(env_dat[occ == 1, , , drop=FALSE],
                               mask = mask,
                               num_starts = 100)

bdry_optim_results <- list()
for (j in 1:nrow(new_starts))
{
  bdry_optim_results[[j]] <- optim(par = new_starts[j,],
                                  fn = loglik_math,
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

    ##      bestlogliks convergences parms_dists      mu1      mu2 sigltil1  sigltil2
    ## [1,]    1174.379            0    0.000000 29.65082 50.55910 6.867986 0.3176915
    ## [2,]    1174.389            0    0.849468 29.57776 49.78176 6.802112 0.3203443
    ## [3,]    1174.404            0    3.430585 29.48900 47.48364 6.662469 0.3241108
    ## [4,]    1174.408            0    2.895631 29.47144 47.91115 6.650329 0.3211054
    ## [5,]    1174.416            0    3.099032 29.46382 47.78718 6.679151 0.3180745
    ## [6,]    1174.425            0    2.639281 29.50804 48.07722 6.617432 0.3225611
    ## [7,]    1174.427            0    3.479017 29.41911 47.38049 6.604492 0.3214007
    ## [8,]    1174.431            0    3.778824 29.41904 47.12596 6.598020 0.3189741
    ##      sigrtil1  sigrtil2      ctil        pd     o_mat1    o_mat2     o_mat3
    ## [1,]      Inf 0.6353020 -13.85134 0.5384179 0.08257498 0.9965849 -0.9965849
    ## [2,]      Inf 0.6332287 -13.51775 0.5389079 0.08254960 0.9965870 -0.9965870
    ## [3,]      Inf 0.6364821 -12.34129 0.5438352 0.08590844 0.9963030 -0.9963030
    ## [4,]      Inf 0.6362675 -12.69390 0.5412129 0.08401464 0.9964645 -0.9964645
    ## [5,]      Inf 0.6456402 -12.47846 0.5432686 0.08353679 0.9965047 -0.9965047
    ## [6,]      Inf 0.6281655 -12.96645 0.5395354 0.08515733 0.9963675 -0.9963675
    ## [7,]      Inf 0.6380317 -12.45674 0.5422012 0.08372054 0.9964893 -0.9964893
    ## [8,]      Inf 0.6435340 -12.28970 0.5429726 0.08401353 0.9964646 -0.9964646
    ##          o_mat4
    ## [1,] 0.08257498
    ## [2,] 0.08254960
    ## [3,] 0.08590844
    ## [4,] 0.08401464
    ## [5,] 0.08353679
    ## [6,] 0.08515733
    ## [7,] 0.08372054
    ## [8,] 0.08401353

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

pnames <-  names(make_mask_names(2))
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
  all_profiles[[counter]] <-  profile_likelihood(
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
likelihood. The `mu1`, `mu2`, and `ctil` profiles, in particular, show
problems. These results suggest a ridge in the likelihood surface that
may rise asymptotically along some path in parameter space for which
`mu2` is increasing. We have not identified a maximum of the likelihood
function, and it looks as though there may not be one if the increase is
indeed asymptotic.

We look at what the growth-environment function looks like for the best
parameters we have found so far, maybe that will give some insight into
what is going wrong:

``` r

param_list <- bdry_optim_results[[1]]$par
param_list["sigrtil1"] <- Inf
param_list <- param_list[names(make_mask_names(2))]
param_list_bio <- math_to_bio(param_list)
interpret_parameters(param_list = param_list_bio,
                          plot_indices = c(1,2), env_dat=env_dat, occ = occ)
```

![](04-unuseable_models_files/figure-html/look_at_inferred_ge_func-1.png)

``` r

param_list_bio$mu
```

    ## [1] 29.65082 50.55910

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
    ## [1,]    1255.413            0 0.0000000000 21.33098 2.399451  761754.5861
    ## [2,]    1255.413            0 0.0005819941 21.33117 2.399351 5079516.6690
    ## [3,]    1255.413            0 0.0007123127 21.33052 2.399294    1840.6994
    ## [4,]    1255.413            0 0.0010142393 21.33094 2.399413     987.1259
    ## [5,]    1255.413            0 0.0012137184 21.33070 2.399331     848.0050
    ## [6,]    1255.413            0 0.0029040052 21.33158 2.399667     352.5439
    ## [7,]    1255.413            0 0.0124616859 21.34342 2.403143 1097898.4436
    ## [8,]    1255.413            0 0.0116930564 21.34101 2.401600     226.9582
    ##           ctil        pd o_mat
    ## [1,] -1.806544 0.6390076     1
    ## [2,] -1.807090 0.6389271     1
    ## [3,] -1.806581 0.6389937     1
    ## [4,] -1.806601 0.6390110     1
    ## [5,] -1.806596 0.6390001     1
    ## [6,] -1.806370 0.6390622     1
    ## [7,] -1.806909 0.6392376     1
    ## [8,] -1.810607 0.6386465     1

One can see the boundary model with `sigrtil` set to `Inf` should be
considered:

``` r

env_dat <- env_array[,,models[i,]==1,drop=FALSE]
mask <- c(sigrtil1 = Inf)
new_starts <- start_parms(env_dat[occ==1,,,drop=FALSE],mask=mask,
                               num_starts=100)

bdry_optim_results5 <- list()
for (j in 1:nrow(new_starts))
{
  bdry_optim_results5[[j]] <- optim(par=new_starts[j,],fn=loglik_math,
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
pnames <- names(make_mask_names(1))
pnames <- pnames[pnames!="sigrtil1"]
linc <- c(0.05,0.025,0.05,0.05)
rinc <- c(0.15,0.025,0.05,0.05)
for (counter in 1:length(pnames))
{
  all_profiles[[counter]] <-  profile_likelihood(
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
