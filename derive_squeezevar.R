true_tau_zero <- 0.6
true_df_zero <- 7.1
n <- 20000
df_g <- 5

# sigma2_g <- extraDistr::rinvchisq(n = n, nu = true_df_zero, tau = true_tau_zero * true_df_zero)
sigma2_g <- 1/(rchisq(n = n, df = true_df_zero))* true_tau_zero * true_df_zero
#sigma2_g <- extraDistr::rinvchisq(n = n, nu = true_df_zero) * true_tau_zero * true_df_zero
s2_g <- rchisq(n = n, df = df_g) * sigma2_g / df_g
hist(s2_g)

# rinvchisqscaled <- function(n, scale, df){
#   # extraDistr::rinvchisq(n, nu =df) / scale
#   1 / rchisq(n, df = df * 2) / scale
# }

dinvchisqscaled <- function(x, scale, df){
  extraDistr::dinvchisq(x, nu = df, tau = (scale / df))
}

dfscale <- function(x, scale, df1, df2, log = FALSE){
  if(log){
    df(x*scale, df1 = df1, df2 = df2, log = TRUE) + log(scale)
  }else{
    df(x*scale, df1 = df1, df2 = df2)*scale
  }
}

hist(sigma2_g, probability = TRUE, breaks = 2500, xlim = c(0, 15))
lines(col = "red", x = seq(0, 15, 0.01), 
      y = dinvchisqscaled(x = seq(0, 15, 0.01), scale = (true_df_zero * true_tau_zero), 
                          df = true_df_zero))


hist(s2_g, probability = TRUE, breaks = 250)
lines(col = "red", x = seq(0, 15, 0.01),
      y = dfscale(x = seq(0, 15, 0.01), 1/true_tau_zero, df_g, true_df_zero))


#chi2 <- 1/rchisq(100000, df = true_df_zero)
sigma2_g <- 1/(rchisq(n = n, df = true_df_zero)) * true_tau_zero  * true_df_zero
#chi2 <- sigma2_g[1:100000]
# chi1 <- rchisq(100000, df = 2)
#s2_g <- rchisq(n = n, df = df_g) * sigma2_g / df_g
s2_g <-  rchisq(n = n, df = df_g) / df_g * sigma2_g 

# f <- chi1*chi2 *(true_df_zero/2)  * true_tau_zero
f <- s2_g#*chi2  / sigma2_g
hist(f, probability = TRUE, breaks = 500, xlim = c(0, 10))
lines(col = "red", x = seq(0, 150, 0.01),
      y = dfscale(x = seq(0, 150, 0.01), 1/true_tau_zero, df_g, true_df_zero))

optim(c(tau0 = 1, df0 = 1.4), function(par){
  tau0 <- par[1]
  df0 <- par[2]
  - sum(dfscale(x = s2_g, 1/tau0, df_g, df0, log = TRUE))
})

limma:::fitFDist(s2_g, df_g)


### Example
df_1 <-  14
df_2 <- 7
n = 500
chi_sq_numerator <- rchisq(n, df = df_1)
chi_sq_denominator <- rchisq(n, df = df_2)

hist(chi_sq_numerator, breaks = 100, probability = TRUE)
xg <- seq(0, 120, by = 0.1)
lines(xg, dchisq(xg, df = df_1), col = "red")

hist(chi_sq_denominator, breaks = 100, probability = TRUE)
lines(xg, dchisq(xg, df = df_2), col = "red")

f_standard <- (chi_sq_numerator/df_1)/(chi_sq_denominator/df_2)

hist(f_standard, breaks = 100, probability = TRUE)
lines(xg, df(xg, df1 = df_1, df2 = df_2), col = "red")

chi_sq_denominator_squeezeVar <- squeezeVar(chi_sq_denominator, df = df_1)
df_0 <- chi_sq_denominator_squeezeVar$df.prior
squeezed_chi_sq_denominator <- chi_sq_denominator_squeezeVar$var.post

f_moderated <- (chi_sq_numerator/df_1)/(squeezed_chi_sq_denominator/(df_2))

hist(f_moderated, breaks = 100, probability = TRUE)
lines(xg, df(xg, df1 = df_1, df2 = df_0 + df_2), col = "red")




