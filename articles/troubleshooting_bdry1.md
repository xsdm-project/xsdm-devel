# Troubleshooting: boundary models, 1

\\ \newcommand{\mean}\[1\]{\overline{#1}} \newcommand{\var}{\text{var}}
\newcommand{\cov}{\text{cov}} \newcommand{\cor}{\text{cor}}
\newcommand{\Rp}{\text{Re}} \newcommand{\E}{\text{E}}
\newcommand{\ltsgr}{\text{ltsgr}} \newcommand{\expit}{\text{expit}}
\newcommand{\logit}{\text{logit}} \\

**Abstract.** We here show a troublshooting example in the case where
the best model considered does not show adequate evidence that the
likelihood function was fully optimized, the best parameters obtained by
optimization show \\p_d\\ close to \\1\\, and the profile for \\p_d\\ is
not dome shaped. This is a common case. The boundary model with \\p_d =
1\\ is used. This example is based on occurrence data from GBIF for
*Blarina carolinensis*, the southern short-tailed shrew.

The southern short-tailed shrew, *Blarina carolinensis*, is found in the
southeastern United States.

We start by loading in the data:

``` r

library(xsdm)
env_array <- example_3$env_array
dim(env_array)
```

    ## [1] 1156   39    6

``` r

dimnames(env_array)[[3]]
```

    ## [1] "BIO01" "BIO10" "BIO11" "BIO12" "BIO16" "BIO17"

``` r

occ <- example_3$occ_vec
length(occ)
```

    ## [1] 1156

Here, there are 6 environmental variables recorded for 39 years in 1156
locations, with accompanying detections and pseudo-absences in the
variable `occ`. The first three environmental variables (BIO1, BIO10,
BIO11) are temperature variables, and the last three (BIO12, BIO16,
BIO17) are precipitation variables. BIO1 is mean annual temperature,
BIO10 is mean temperature of the warmest quarter, BIO11 is mean
temperature of the coldest quarter, BIO12 is annual precipitation, BIO16
is precipitation of the wettest quarter, and BIO17 is precipitation of
the driest quarter.

Now look at the distributions of values of environmental variables to
make sure they are not on wildly different scales, which would cause
problems for optimization:

``` r

apply(FUN=quantile, X=env_array, MARGIN=3,prob=c(.025,.25,.5,.75,.975))
```

    ##          BIO01    BIO10     BIO11     BIO12    BIO16     BIO17
    ## 2.5%  11.61870 21.70901  1.477658  76.36044  9.12790  2.133369
    ## 25%   15.31648 25.24791  6.289062 108.62843 13.18716  4.505094
    ## 50%   17.16010 26.45451  9.055910 126.29746 15.76004  5.947245
    ## 75%   19.20414 27.43751 11.991494 146.79118 18.83666  7.377723
    ## 97.5% 22.55861 29.21995 17.613174 190.28207 27.04961 10.422263

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
  env_dat <- env_array[,,models[i,]==1,drop=FALSE]
  starts <- xsdm::start_parms(env_dat,num_starts=25)
  all_optim_results <- list()
  for (j in 1:nrow(starts))
  {
    all_optim_results[[j]] <- optim(par=starts[j,],fn=xsdm::loglik_math,
                                   method="BFGS",
                                   env_dat=env_dat, occ=occ,negative=TRUE,
                                   control=list(trace=0))
  }
  all_model_results[[i]] <- all_optim_results
}
```

Rank the models by BIC, bearing in mind that we’ve been working with the
negative of the likelihood:

``` r

model_BICs <- sapply(X=all_model_results,
                      FUN=function(x){
                        best_loglik = min(sapply(X=x, FUN=function(y){y$value}))
                        num_parms = length(x[[1]]$par)
                        n = length(occ)
                        BIC = 2*best_loglik + num_parms*log(n)
                        return(BIC)
                      }
                    )
```

Also by AIC:

``` r

model_AICs <- sapply(X=all_model_results,
                      FUN=function(x){
                        best_loglik = min(sapply(X=x, FUN=function(y){y$value}))
                        num_parms = length(x[[1]]$par)
                        AIC = 2*best_loglik + 2*num_parms
                        return(AIC)
                      }
                    )
inds <- order(model_BICs)
rbind(model_BICs[inds],model_AICs[inds])
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
    ## [1,] 1165.031 1172.509 1175.884 1190.920 1193.879 1194.610 1196.464 1201.064
    ## [2,] 1119.557 1127.035 1150.620 1145.446 1148.404 1149.136 1171.200 1155.589
    ##          [,9]    [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 1201.254 1202.077 1210.059 1216.145 1288.876 1292.384 1310.530
    ## [2,] 1155.779 1156.602 1164.584 1190.882 1263.613 1267.121 1285.266

``` r

plot(model_BICs,model_AICs,type="p",xlab="BIC",ylab="AIC")
```

![](troubleshooting_bdry1_files/figure-html/rank_by_AIC-1.png)

``` r

order(model_BICs)
```

    ##  [1]  9 12  2  7  8 15  1 11 10 13 14  3  5  4  6

``` r

order(model_AICs)
```

    ##  [1]  9 12  7  8 15  2 11 10 13 14  1  3  5  4  6

Looks like in this case the AIC and the BIC are pretty aligned with each
other. Look at the best two model models:

``` r

models[order(model_BICs)[1:2],]
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]    1    0    0    0    0    1
    ## [2,]    0    1    0    0    0    1

So these are two-variable models.

Let’s use the best model. Optimize it a bit harder:

``` r

i <- 9
env_dat <- env_array[,,models[i,]==1,drop=FALSE]
starts <- xsdm::start_parms(env_dat,num_starts=100)
best_model_results <- list()
for (j in 1:nrow(starts))
{
  best_model_results[[j]] <- optim(par=starts[j,],fn=xsdm::loglik_math,
                                 method="BFGS",
                                 env_dat=env_dat, occ=occ,negative=TRUE,
                                 control=list(trace=0))
}
values <- sapply(X=best_model_results, FUN=function(y){y$value})
inds <- order(values)
best_model_results <- best_model_results[inds]
min(sapply(X=all_model_results[[i]], FUN=function(y){y$value}))
```

    ## [1] 550.7784

``` r

best_model_results[[1]]$value
```

    ## [1] 550.7784

Pretty similar, so we HAD pretty much fully optimized before.

Now have a look at the result for this model.

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

  #look at the best 5 results in parameter space
  bestparms <- sapply(X=parms_dists_to_best, FUN=function(x){unlist(x$representative)})

  #put it all together
  return(rbind(bestlogliks,convergences,parms_dists,bestparms))
}

h <- examine_optim_results(best_model_results)
t(h[,1:8])
```

    ##      bestlogliks convergences  parms_dists      mu1      mu2 sigltil1  sigltil2
    ## [1,]    550.7784            0 0.000000e+00 14.26892 3.549587 4.602682 0.4576370
    ## [2,]    550.7784            0 3.461061e-05 14.26890 3.549594 4.602665 0.4576380
    ## [3,]    550.7784            0 1.699042e-04 14.26883 3.549640 4.602578 0.4576106
    ## [4,]    550.7784            0 2.487303e-04 14.26900 3.549521 4.602938 0.4576604
    ## [5,]    550.7784            0 4.220242e-05 14.26893 3.549586 4.602695 0.4576419
    ## [6,]    550.7784            0 1.420623e-04 14.26896 3.549556 4.602595 0.4576442
    ## [7,]    550.7784            0 4.031916e-05 14.26892 3.549582 4.602692 0.4576393
    ## [8,]    550.7784            0 2.497969e-05 14.26892 3.549588 4.602676 0.4576396
    ##       sigrtil1 sigrtil2       ctil        pd     o_mat1     o_mat2    o_mat3
    ## [1,] 0.1481576 4.892043 -0.1921917 0.9999952 -0.2429301 -0.9700438 0.9700438
    ## [2,] 0.1481577 4.892054 -0.1921622 0.9999939 -0.2429355 -0.9700424 0.9700424
    ## [3,] 0.1481588 4.892075 -0.1922056 0.9999936 -0.2429382 -0.9700418 0.9700418
    ## [4,] 0.1481535 4.892037 -0.1921324 0.9999912 -0.2429303 -0.9700438 0.9700438
    ## [5,] 0.1481571 4.891990 -0.1921704 0.9999907 -0.2429351 -0.9700426 0.9700426
    ## [6,] 0.1481548 4.892035 -0.1922131 0.9999907 -0.2429281 -0.9700443 0.9700443
    ## [7,] 0.1481568 4.892010 -0.1921813 0.9999900 -0.2429346 -0.9700427 0.9700427
    ## [8,] 0.1481574 4.892031 -0.1921739 0.9999881 -0.2429351 -0.9700425 0.9700425
    ##          o_mat4
    ## [1,] -0.2429301
    ## [2,] -0.2429355
    ## [3,] -0.2429382
    ## [4,] -0.2429303
    ## [5,] -0.2429351
    ## [6,] -0.2429281
    ## [7,] -0.2429346
    ## [8,] -0.2429351

This looks good, in the sense that it looks like multiple optimizations
arrived at about the same place in parameter space, which is some
evidence that we may have found the global maximum of the likelihood
function The only problem is that \\p_d\\ is very close to \\1\\.

Let’s profile this model and see what we get:

``` r

pnames <- names(xsdm::make_mask_names(2))

all_profiles <- list()
linc <- c(rep(0.05,8),0.01)
rinc <- c(rep(0.05,8),0.01)
for (counter in 1:9)
{
  all_profiles[[counter]] <- xsdm::profile_likelihood(
                              profile_parameter=pnames[counter],
                              increment_left=linc[counter],
                              increment_right=rinc[counter],
                              num_steps_left=50,
                              num_steps_right=50,
                              alpha=0.95,
                              optim_param_vector=best_model_results[[1]]$par,
                              env_dat=env_dat,
                              occ=occ,
                              mask=NULL,
                              num_threads=6
                            )
}
names(all_profiles) <- pnames
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
plot_tool(all_profiles,1)
plot_tool(all_profiles,2)
plot_tool(all_profiles,3)
plot_tool(all_profiles,4)
plot_tool(all_profiles,5)
plot_tool(all_profiles,6)
plot_tool(all_profiles,7)
plot_tool(all_profiles,8)
plot_tool(all_profiles,9)
```

![](troubleshooting_bdry1_files/figure-html/plot_profiles_first-1.png)

So yes, there is the expected problem with pd. This is a common problem
and it manifests in the manner illustrated here.

The solution is to use the boundary model with \\p_d = 1\\:

``` r

env_dat <- env_array[,,models[i,]==1,drop=FALSE]
dim(env_dat)
```

    ## [1] 1156   39    2

``` r

mask <- c(pd=Inf) #use Inf because masks are given on the math scale
mask
```

    ##  pd 
    ## Inf

``` r

new_starts <- xsdm::start_parms(env_dat[occ==1,,,drop=FALSE],
                                  mask=mask,num_starts=100)
head(new_starts)
```

    ## # A tibble: 6 × 8
    ##     mu1   mu2 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil o_par1
    ##   <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
    ## 1  15.6  5.66   0.232    1.09      0.562    1.27  -1.19    3.68
    ## 2  18.1  8.30   0.925    0.397     1.26     0.572 -0.618  -5.74
    ## 3  19.3  4.34   0.579    0.744     1.60     0.919 -0.905   8.39
    ## 4  16.9  6.98  -0.115    0.0504    0.908    0.226 -1.48   -1.03
    ## 5  17.5  3.68   0.752    1.26      1.43     0.399 -1.34   -8.10
    ## 6  19.9  6.32   0.0587   0.570     0.735    1.09  -0.762   1.33

``` r

bdry_optim_results <- list()
for (j in 1:nrow(new_starts))
{
  bdry_optim_results[[j]] <- optim(par=new_starts[j,],fn=xsdm::loglik_math,
                                method="BFGS",
                                env_dat=env_dat,occ=occ,mask=mask,negative=TRUE,
                                control=list(trace=0,maxit=500))
}
```

Have a look at these results:

``` r

h <- examine_optim_results(bdry_optim_results,mask=mask)
t(h[,1:10])
```

    ##       bestlogliks convergences  parms_dists      mu1      mu2 sigltil1
    ##  [1,]    550.7784            0 0.000000e+00 14.26893 3.549582 4.602702
    ##  [2,]    550.7784            0 1.524294e-04 14.26900 3.549559 4.602727
    ##  [3,]    550.7784            0 3.314862e-04 14.26895 3.549538 4.602675
    ##  [4,]    550.7784            0 9.525015e-05 14.26898 3.549575 4.602621
    ##  [5,]    550.7784            0 2.851435e-04 14.26895 3.549554 4.602495
    ##  [6,]    550.7784            0 1.052428e-04 14.26895 3.549580 4.602604
    ##  [7,]    550.7784            0 1.609023e-04 14.26900 3.549560 4.602611
    ##  [8,]    550.7784            0 2.161842e-04 14.26893 3.549569 4.602639
    ##  [9,]    550.7784            0 9.790124e-05 14.26892 3.549588 4.602698
    ## [10,]    550.7784            0 1.415520e-04 14.26892 3.549583 4.602651
    ##        sigltil2  sigrtil1 sigrtil2       ctil pd     o_mat1     o_mat2
    ##  [1,] 0.4576392 0.1481599 4.892063 -0.1921664  1 -0.2429293 -0.9700440
    ##  [2,] 0.4576672 0.1481599 4.892012 -0.1921522  1 -0.2429288 -0.9700441
    ##  [3,] 0.4576516 0.1481528 4.892021 -0.1921763  1 -0.2429307 -0.9700436
    ##  [4,] 0.4576555 0.1481604 4.892036 -0.1921767  1 -0.2429312 -0.9700435
    ##  [5,] 0.4576450 0.1481538 4.892029 -0.1922115  1 -0.2429313 -0.9700435
    ##  [6,] 0.4576506 0.1481580 4.892035 -0.1921768  1 -0.2429314 -0.9700435
    ##  [7,] 0.4576662 0.1481612 4.892100 -0.1921539  1 -0.2429290 -0.9700441
    ##  [8,] 0.4576433 0.1481552 4.892042 -0.1921706  1 -0.2429320 -0.9700433
    ##  [9,] 0.4576423 0.1481578 4.892051 -0.1921685  1 -0.2429328 -0.9700431
    ## [10,] 0.4576357 0.1481569 4.892026 -0.1921776  1 -0.2429298 -0.9700439
    ##          o_mat3     o_mat4
    ##  [1,] 0.9700440 -0.2429293
    ##  [2,] 0.9700441 -0.2429288
    ##  [3,] 0.9700436 -0.2429307
    ##  [4,] 0.9700435 -0.2429312
    ##  [5,] 0.9700435 -0.2429313
    ##  [6,] 0.9700435 -0.2429314
    ##  [7,] 0.9700441 -0.2429290
    ##  [8,] 0.9700433 -0.2429320
    ##  [9,] 0.9700431 -0.2429328
    ## [10,] 0.9700439 -0.2429298

This looks like enough of them converged to the same thing to say that
it looks like I successfully optimized. Get the likelihood:

``` r

values <- sapply(X=bdry_optim_results, FUN=function(y){y$value})
inds <- order(values)
bdry_optim_results <- bdry_optim_results[inds]
best_model_results[[1]]$value
```

    ## [1] 550.7784

``` r

bdry_optim_results[[1]]$value
```

    ## [1] 550.7784

So the likelihood is the same as for the previous (non-boundary) model.
But we have one less parameter that has been fitted. Get the AIC and BIC
to see the effects:

``` r

AIC <- unname(2*bdry_optim_results[[1]]$value+
               2*length(bdry_optim_results[[1]]$par))
AIC_old <- unname(2*best_model_results[[1]]$value+
                   2*length(best_model_results[[1]]$par))
AIC
```

    ## [1] 1117.557

``` r

AIC_old
```

    ## [1] 1119.557

``` r

AIC-AIC_old
```

    ## [1] -2.000041

``` r

BIC <- unname(2*bdry_optim_results[[1]]$value+
               log(length(occ))*length(bdry_optim_results[[1]]$par))
BIC_old <- unname(2*best_model_results[[1]]$value+
               log(length(occ))*length(best_model_results[[1]]$par))
BIC
```

    ## [1] 1157.979

``` r

BIC_old
```

    ## [1] 1165.031

``` r

BIC-BIC_old
```

    ## [1] -7.052762

``` r

log(length(occ))
```

    ## [1] 7.052721

So the AIC and BIC are better for the boundary model compared to the
earlier model, and by the expected amounts. We go with the simpler
model.

Let’s profile this boundary model and see what we get:

``` r

pnames <- names(xsdm::make_mask_names(2))
pnames <- pnames[pnames!="pd"]
mask
```

    ##  pd 
    ## Inf

``` r

all_bdry_profiles <- list()
linc <- c(rep(0.05,7),0.01)
rinc <- c(rep(0.05,7),0.01)
for (counter in 1:8)
{
  all_bdry_profiles[[counter]] <- xsdm::profile_likelihood(
                              profile_parameter=pnames[counter],
                              increment_left=linc[counter],
                              increment_right=rinc[counter],
                              num_steps_left=50,
                              num_steps_right=50,
                              alpha=0.95,
                              optim_param_vector=bdry_optim_results[[1]]$par,
                              env_dat=env_dat,
                              occ=occ,
                              mask=mask,
                              num_threads=6
                            )
}
names(all_bdry_profiles) <- pnames
```

Now plot these profiles:

``` r

par(mfrow=c(3,3))
plot_tool(all_bdry_profiles,1)
plot_tool(all_bdry_profiles,2)
plot_tool(all_bdry_profiles,3)
plot_tool(all_bdry_profiles,4)
plot_tool(all_bdry_profiles,5)
plot_tool(all_bdry_profiles,6)
plot_tool(all_bdry_profiles,7)
plot_tool(all_bdry_profiles,8)
```

![](troubleshooting_bdry1_files/figure-html/plot_profiles-1.png)

These are dome-shaped, suggesting the likelihood function is
well-behaved in the vacinity of the maximum we found, and inferences can
be made. Note these profiles look essentially the same as those obtained
previously, except \\\vec{\sigma}\_L\\ and \\\vec{\sigma}\_R\\ have been
switched (which can happen, given redundancy in parameters explained in
“How to fit xsdm models with species occurrence data using xsdm”). The
overall lesson is, when optimizing the xsdm likelihood function gives a
best model with \\p_d\\ close to \\1\\ and a non-dome-shaped profile for
\\p_d\\, use the boundary model with \\p_d = 1\\.
