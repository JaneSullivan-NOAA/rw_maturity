library(dplyr)
library(TMB)
library(reshape2) # for melt()
library(lattice) # for matrix plotting
library(ggplot2)
precompile() # precompiling of tmb libraries, slow the first time, faster compilation after that

# dyn.unload(dynlib("maturity_rw"))
compile("rw_maturity.cpp")
dyn.load(dynlib("rw_maturity"))

#  ----
set.seed(123)

# simulation values for logistic maturity curve with random walk on a50 (median, age at 50%
# maturity) and kmat (slope)
a50_init = 9
a50_sd = 0.2
kmat_init = 0.8
kmat_sd = 0.05
nage = 20
nyear = 100

# random walk
a50_init_vec <- a50_init + cumsum(rnorm(n = nyear, mean=0, sd = a50_sd))
kmat_init_vec <- kmat_init + cumsum(rnorm(n = nyear, mean=0, sd = kmat_sd))
a50_init_vec
kmat_init_vec
ages = t(matrix(1:nage, nage, nyear))
mat <- 1 / (1 + exp(-kmat_init_vec * (ages - a50_init_vec)))

lattice::levelplot(mat, xlab='Year index', ylab='Age')
plot(1:nyear, kmat_init_vec)
plot(1:nyear, a50_init_vec)

# simulate data (N=total fish, Y=mature fish)
N <- matrix(sample(50:100, length(ages), replace=TRUE), nyear, nage)
Y <- matrix(rbinom(length(N), N, mat), nyear, nage)

# create parameter, data, and random objects
input = list(data=list(), par=list())
input$par$log_a50_sd = log(a50_sd)
input$par$log_kmat_sd = log(kmat_sd)
input$par$log_a50_vec = rep(log(a50_init), nyear)
input$par$log_kmat_vec = rep(log(kmat_init), nyear)

random = c("log_a50_vec", "log_kmat_vec")

input$data$Y = Y
input$data$N = N
input$data$age_obs = ages
input$data$max_age = max(ages)
map = NULL

mod = TMB::MakeADFun(input$data, input$par,
                     map = map, random = random,
                     DLL = "rw_maturity")

opt = nlminb(mod$par, mod$fn, mod$gr)
mod$rep = mod$report()
mod$sdrep = sdreport(mod)
# summary(mod$sdrep)
mod$sdrep # fixed effects of the process error estimates for a50 and kmat

lattice::levelplot(mod$rep$mat, xlab='Year index', ylab='Age')
lattice::levelplot(mat, xlab='Year index', ylab='Age')
plot(1:nyear, kmat_init_vec)
points(1:nyear, as.list(mod$sdrep, "Estimate", report = TRUE)$kmat_vec, col = 'red', 'l')
plot(1:nyear, a50_init_vec)
points(1:nyear, as.list(mod$sdrep, "Estimate", report = TRUE)$a50_vec, col = 'red', 'l')

logit_mat = 

annual_matdf <- mod$rep$mat %>% 
  reshape2::melt() %>% 
  setNames(c('year', 'age', 'pmat')) 

matdf <- data.frame(age = 1:nage, 
                    pmat = 1/(1+exp(-as.list(mod$sdrep, "Estimate", report = TRUE)$logit_mean_mat)))

ggplot() +
  geom_line(data = annual_matdf, aes(x = age, y = pmat, col = year, group = year)) +
  geom_line(data = matdf, aes(x = age, y = pmat), size = 3) 
