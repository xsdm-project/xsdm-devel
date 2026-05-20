# Boundary models

\\ \newcommand{\mean}\[1\]{\overline{#1}} \newcommand{\var}{\text{var}}
\newcommand{\cov}{\text{cov}} \newcommand{\cor}{\text{cor}}
\newcommand{\Rp}{\text{Re}} \newcommand{\E}{\text{E}}
\newcommand{\ltsgr}{\text{ltsgr}} \newcommand{\expit}{\text{expit}}
\newcommand{\logit}{\text{logit}} \\

**Abstract.** We here show an example of the not uncommon case where the
seemingly best initial model considered does not show adequate evidence
that the likelihood function was fully optimized; the best parameters
obtained by optimization show \\p_d\\ close to \\1\\; and the profile
for \\p_d\\ is not dome shaped. The boundary model with \\p_d = 1\\ can
then be used. This example is based on occurrence data from GBIF for
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
make sure they are not on wildly different scales, which could cause
problems for optimization:

``` r

apply(FUN=quantile, X=env_array, MARGIN=3, prob=c(.025,.25,.5,.75,.975))
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
  starts <- start_parms(env_dat[occ==1,,,drop=FALSE],num_starts=25)
  all_optim_results <- list()
  for (j in 1:nrow(starts))
  {
    all_optim_results[[j]] <- optim(par=starts[j,],fn=loglik_math,
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

Also by AIC, and compare the two:

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
    ## [1,] 1165.031 1172.509 1175.884 1190.924 1193.879 1194.610 1196.464 1201.064
    ## [2,] 1119.557 1127.035 1150.620 1145.450 1148.404 1149.136 1171.200 1155.590
    ##          [,9]    [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 1201.256 1202.077 1210.059 1216.145 1288.876 1292.384 1310.530
    ## [2,] 1155.781 1156.602 1164.584 1190.881 1263.613 1267.121 1285.266

``` r

plot(model_BICs,model_AICs,type="p",xlab="BIC",ylab="AIC")
```

![](03-boundary_models_files/figure-html/rank_by_AIC-1.png)

``` r

order(model_BICs)
```

    ##  [1]  9 12  2  7  8 15  1 11 10 13 14  3  5  4  6

``` r

order(model_AICs)
```

    ##  [1]  9 12  7  8 15  2 11 10 13 14  1  3  5  4  6

Looks like in this case the AIC and the BIC are reasonably aligned with
each other. Look at the best two models:

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
starts <- start_parms(env_dat[occ==1,,,drop=FALSE],num_starts=100)
best_model_results <- list()
for (j in 1:nrow(starts))
{
  best_model_results[[j]] <- optim(par=starts[j,],fn=loglik_math,
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
      dist_between_params(
        x$par,
        best_parms_math,
        mask=mask,
        give_closest_rep=TRUE)
    }
  )
  parms_dists <- sapply(X=parms_dists_to_best, FUN=function(x){x$distance})

  #get the parameters from each optimization run which are closest to the first set
  bestparms <- sapply(X=parms_dists_to_best, FUN=function(x){unlist(x$representative)})

  #put it all together
  return(rbind(bestlogliks,convergences,parms_dists,bestparms))
}

h <- examine_optim_results(best_model_results)
t(h[,1:8])
```

    ##      bestlogliks convergences  parms_dists      mu1      mu2  sigltil1
    ## [1,]    550.7784            0 0.0000000000 14.26891 3.549614 0.4576326
    ## [2,]    550.7784            0 0.0001105175 14.26894 3.549601 0.4576453
    ## [3,]    550.7784            0 0.0002127497 14.26889 3.549585 0.4576303
    ## [4,]    550.7784            0 0.0001552695 14.26887 3.549602 0.4576217
    ## [5,]    550.7784            0 0.0001600068 14.26892 3.549588 0.4576400
    ## [6,]    550.7784            0 0.0001473951 14.26896 3.549614 0.4576482
    ## [7,]    550.7784            0 0.0001287526 14.26894 3.549578 0.4576435
    ## [8,]    550.7784            0 0.0002180582 14.26896 3.549564 0.4576563
    ##       sigltil2 sigrtil1 sigrtil2       ctil        pd    o_mat1     o_mat2
    ## [1,] 0.1481606 4.891978 4.602757 -0.1921649 0.9999959 0.9700422 -0.2429363
    ## [2,] 0.1481589 4.892003 4.602559 -0.1921996 0.9999953 0.9700424 -0.2429356
    ## [3,] 0.1481560 4.892077 4.602709 -0.1921684 0.9999950 0.9700439 -0.2429296
    ## [4,] 0.1481575 4.892052 4.602718 -0.1921667 0.9999926 0.9700423 -0.2429360
    ## [5,] 0.1481573 4.892028 4.602664 -0.1921785 0.9999904 0.9700426 -0.2429350
    ## [6,] 0.1481631 4.891955 4.602518 -0.1921786 0.9999881 0.9700427 -0.2429343
    ## [7,] 0.1481583 4.892035 4.602685 -0.1921474 0.9999878 0.9700427 -0.2429343
    ## [8,] 0.1481569 4.892021 4.602692 -0.1921702 0.9999839 0.9700430 -0.2429331
    ##         o_mat3    o_mat4
    ## [1,] 0.2429363 0.9700422
    ## [2,] 0.2429356 0.9700424
    ## [3,] 0.2429296 0.9700439
    ## [4,] 0.2429360 0.9700423
    ## [5,] 0.2429350 0.9700426
    ## [6,] 0.2429343 0.9700427
    ## [7,] 0.2429343 0.9700427
    ## [8,] 0.2429331 0.9700430

This looks good, in the sense that it looks like multiple optimizations
arrived at the same place in parameter space, which is some evidence
that we may have found the global maximum of the likelihood function.
The only problem is that \\p_d\\ is very close to \\1\\.

Let’s profile this model and see what we get:

``` r

pnames <- names(make_mask_names(2))

all_profiles <- list()
linc <- c(rep(0.05,8),0.01)
rinc <- c(rep(0.05,8),0.01)
for (counter in 1:9)
{
  all_profiles[[counter]] <- profile_likelihood(
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

![](03-boundary_models_files/figure-html/plot_profiles_first-1.png)

So yes, the expected problem with `pd` is appearing. We cannot use a
model for which profiles do not show a dome-like pattern. This is a
common problem and it manifests in the manner illustrated here.

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

new_starts <- start_parms(env_dat[occ==1,,,drop=FALSE],
                                  mask=mask,num_starts=100)
head(new_starts)
```

    ## # A tibble: 6 × 8
    ##     mu1   mu2 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil o_par1
    ##   <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
    ## 1  15.6  7.89   0.297   1.05       1.19     0.702 -1.32   -3.39
    ## 2  18.1  5.25   0.990   0.354      0.497    1.40  -0.744   6.04
    ## 3  19.3  6.57  -0.0496  0.700      1.54     0.356 -1.61   -8.10
    ## 4  16.9  3.92   0.644   0.00711    0.843    1.05  -1.03    1.33
    ## 5  17.5  5.91   0.124   0.874      1.36     0.182 -1.46   -5.74
    ## 6  19.9  8.55   0.817   0.180      0.670    0.875 -0.887   3.68

``` r

bdry_optim_results <- list()
for (j in 1:nrow(new_starts))
{
  bdry_optim_results[[j]] <- optim(par=new_starts[j,],fn=loglik_math,
                                method="BFGS",
                                env_dat=env_dat,occ=occ,mask=mask,negative=TRUE,
                                control=list(trace=0,maxit=500))
}
```

Have a look at these results:

``` r

h <- examine_optim_results(bdry_optim_results,mask=mask)
t(h[,1:8])
```

    ##      bestlogliks convergences  parms_dists      mu1      mu2  sigltil1 sigltil2
    ## [1,]    550.7784            0 0.000000e+00 14.26898 3.549565 0.1481579 4.892005
    ## [2,]    550.7784            0 4.800380e-05 14.26896 3.549569 0.1481577 4.892013
    ## [3,]    550.7784            0 2.621185e-05 14.26898 3.549556 0.1481583 4.892046
    ## [4,]    550.7784            0 1.251270e-04 14.26895 3.549562 0.1481553 4.892064
    ## [5,]    550.7784            0 1.242072e-04 14.26905 3.549547 0.1481583 4.891977
    ## [6,]    550.7784            0 6.697159e-05 14.26896 3.549563 0.1481568 4.892015
    ## [7,]    550.7784            0 1.334626e-04 14.26898 3.549544 0.1481556 4.892069
    ## [8,]    550.7784            0 1.228045e-04 14.26892 3.549584 0.1481572 4.892030
    ##      sigrtil1  sigrtil2       ctil pd    o_mat1    o_mat2     o_mat3    o_mat4
    ## [1,] 4.602648 0.4576596 -0.1921786  1 0.2429296 0.9700439 -0.9700439 0.2429296
    ## [2,] 4.602687 0.4576510 -0.1921726  1 0.2429303 0.9700437 -0.9700437 0.2429303
    ## [3,] 4.602512 0.4576606 -0.1921872  1 0.2429291 0.9700440 -0.9700440 0.2429291
    ## [4,] 4.602624 0.4576539 -0.1921786  1 0.2429291 0.9700440 -0.9700440 0.2429291
    ## [5,] 4.602488 0.4576790 -0.1922139  1 0.2429287 0.9700441 -0.9700441 0.2429287
    ## [6,] 4.602650 0.4576518 -0.1921791  1 0.2429327 0.9700431 -0.9700431 0.2429327
    ## [7,] 4.602353 0.4576443 -0.1922112  1 0.2429266 0.9700447 -0.9700447 0.2429266
    ## [8,] 4.602648 0.4576392 -0.1921772  1 0.2429327 0.9700431 -0.9700431 0.2429327

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

    ## [1] -2.000035

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

    ## [1] -7.052756

``` r

log(length(occ))
```

    ## [1] 7.052721

So the AIC and BIC are better for the boundary model compared to the
earlier model, and by the expected amounts. We go with the boundary
model.

Let’s profile this boundary model and see what we get:

``` r

pnames <- names(make_mask_names(2))
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
  all_bdry_profiles[[counter]] <- profile_likelihood(
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

![](03-boundary_models_files/figure-html/plot_profiles-1.png)

These are sufficiently dome-shaped, suggesting the likelihood function
is sufficiently well-behaved in the vicinity of the maximum we found,
and inferences can be made. Note these profiles look essentially the
same as those obtained previously, except the \\\sigma\\ profiles have
been switched around (which can happen, given redundancy in parameters
explained in “How to fit xsdm models with species occurrence data using
xsdm”). The overall lesson is, when optimizing the xsdm likelihood
function gives a best model with \\p_d\\ close to \\1\\ and a
non-dome-shaped profile for \\p_d\\, use the boundary model with \\p_d =
1\\. We also saw a demonstration of the `xsdm` functionality for fitting
boundary models. Other boundary models can also be fitted. See the
document “Unusable models: when a model does not have a maximum
likelihood” for more examples of fitting boundary models.
