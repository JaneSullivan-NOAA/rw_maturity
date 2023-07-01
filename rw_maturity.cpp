// logistic maturity-at-age model with a random walk on a50 and the slope
// parameter (kmat)
// jane.sullivan@noaa.gov
// june 2023

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Y); //number mature
  DATA_MATRIX(N); //number of mature + not mature
  DATA_MATRIX(age_obs); //age for each observation
  DATA_INTEGER(max_age); //to generate a maturity ogive
  
  PARAMETER_VECTOR(log_a50_vec); // annual a50 value
  PARAMETER_VECTOR(log_kmat_vec); // annual kmat value
  PARAMETER(log_a50_sd); // random walk SD for a50
  PARAMETER(log_kmat_sd); // random walk SD for slope
  
  // set up model dimensions and predicted values
  int nyear = Y.rows();
  int nage = Y.cols();
  Type nll = 0;
  
  vector<Type> a50_vec(nyear);
  vector<Type> kmat_vec(nyear);
  a50_vec.setZero();
  kmat_vec.setZero();
  
  for(int i = 0; i < nyear; i++) {
    a50_vec(i) = exp(log_a50_vec(i));
    kmat_vec(i) = exp(log_kmat_vec(i));
  }
  
  matrix<Type> logit_mat(nyear, nage);
  matrix<Type> mat(nyear, nage);
  vector<Type> logit_mean_mat(max_age);
  vector<Type> mean_mat(max_age);
  
  for(int i = 0; i < nyear; i++) {
    for(int j = 0; j < nage; j++) {
      logit_mat(i,j) = kmat_vec(i) * (age_obs(i,j) - a50_vec(i));
      mat(i,j) = Type(1) / (Type(1) + exp(-logit_mat(i,j)));
      nll -= dbinom(Y(i,j), N(i,j), mat(i,j), 1);   // negative log-likelihood. the 1 means log = TRUE
    }
  }
  
  // predicted mean maturity-at-age (across all years)
  Type a50 = a50_vec.sum() / nyear;
  Type kmat = kmat_vec.sum() / nyear;
  
  for(int i = 0; i < max_age; i++) {
    logit_mean_mat(i) = kmat * (Type(i+1) - a50);
    mean_mat(i) = Type(1) / (Type(1) + exp(-logit_mean_mat(i)));
  }
  
  // random walk likelihood component
  for(int i = 1; i < nyear; i++) {
    nll -= dnorm(log_a50_vec(i-1), log_a50_vec(i), exp(log_a50_sd), 1);
    nll -= dnorm(log_kmat_vec(i-1), log_kmat_vec(i), exp(log_kmat_sd), 1);
  }
  
  REPORT(nll);
  ADREPORT(a50);
  ADREPORT(kmat);
  REPORT(mat);
  ADREPORT(logit_mean_mat);
  ADREPORT(a50_vec);
  ADREPORT(kmat_vec);
  return nll;
}
