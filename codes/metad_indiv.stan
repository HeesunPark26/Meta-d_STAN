// Version 2
// Bayesian estimation of meta-d (not hierarchical)
 // Sampling from Truncated Normal (https://sites.google.com/site/michaelchughessite/research/sampling-from-truncated-normal)
functions {
  real sample_from_tnormal_lb(real a, real mu, real sigma, real u_pr) {
    real phi_a_bar = Phi_approx((a - mu) / sigma);
    real u = phi_a_bar + (1 - phi_a_bar) * u_pr;
    real x = sigma * inv_Phi(u) + mu;
    return x;
  }
  real sample_from_tnormal_ub(real a, real mu, real sigma, real u_pr) {
    real phi_a_bar = Phi_approx((a - mu) / sigma);
    real u = phi_a_bar * u_pr;
    real x = sigma * inv_Phi(u) + mu;
    return x;
  }
  vector normalize(vector x) {
    return x / sum(x);
  }
}

data {
  int<lower=1> N;             //Number of Subjects
  
  // Type 1 data
  int<lower=0> H[N];          // Hit trials (resp=S2, stim=S2)
  int<lower=0> FA[N];         // False Alarm (resp=S2, stim=S1)
  int<lower=0> CR[N];         // Correct rejection(resp=S1, stim=S1)
  int<lower=0> M[N];          // Miss (resp=S1, stim = S2)
  
  //Type 2 data
  int<lower=1> nratings;      //the number of ratings
  int counts[N, 4*nratings];  //count vectors for all subjects, counts = [nR_S1, nR_S2]
  
  // Type 1 SDT parameters
  vector[N] d1;               // Type 1 discriminability
  vector[N] c1;               // Type 1 bias
}

parameters {  
  // Type 2 SDT parameters
  vector[N] meta_d;
  matrix<lower = 0, upper = 1>[N, nratings-1] cS1_p;
  matrix<lower = 0, upper = 1>[N, nratings-1] cS2_p;
}

transformed parameters{
  vector[N] S2mu;                     //mean of S2 (Signal Distribution)
  vector[N] S1mu;                     //mean of S1 (Noise Distribution)
  matrix[N, nratings-1] cS1_raw;      //c2 vectors for S1 (cS1_p -> cS1_raw -> cS1)
  matrix[N, nratings-1] cS2_raw;      //c2 vectors for S2 (cS2_p -> cS2_raw -> cS2)
  matrix[N, nratings-1] cS1;          //c2 vectors (cS1_p -> cS1_raw -> cS1)
  matrix[N, nratings-1] cS2;          //c2 vectors (cS2_p -> cS2_raw -> cS2)
  
  matrix[N, 4] log_lik;

  for (i in 1:N){
    for (j in 1:(nratings-1)) {
      cS1_raw[i, j] = sample_from_tnormal_ub(c1[i], 0., 2., cS1_p[i, j]); //sampling from truncated normal distribution (upper bound)
      cS2_raw[i, j] = sample_from_tnormal_lb(c1[i], 0., 2., cS2_p[i, j]); //sampling from truncated normal distribution (lower bound)
      }
    // specify ordered prior on criteria
    cS1[i, 1:nratings-1] = to_row_vector(sort_asc(cS1_raw[i,]));
    cS2[i, 1:nratings-1] = to_row_vector(sort_asc(cS2_raw[i,]));
  }
  
  // Means of SDT distributions
  S2mu = meta_d/2; 
  S1mu = -meta_d/2;
  

  for (i in 1:N){
    real C_area_rS1;             // Correct rejection
    real I_area_rS1;             // Miss
    real C_area_rS2;             // Hit
    real I_area_rS2;             // False Alarm
    vector[nratings*4] pr;       //vectors for probabilities 
    vector[nratings*4] prT;      //transformed probabilities
        
    // Calculate normalization constants
    C_area_rS1 = Phi_approx(c1[i] - S1mu[i]); 
    I_area_rS1 = Phi_approx(c1[i] - S2mu[i]); 
    I_area_rS2 = 1-C_area_rS1;
    C_area_rS2 = 1-I_area_rS1;

    // Get nC_rS1 probs (Correct Rejection)
    pr[1] = Phi_approx(cS1[i,1] - S1mu[i])/C_area_rS1;  
    for (k in 1:(nratings-2)) {
      pr[k+1] = (Phi_approx(cS1[i, k+1] - S1mu[i])-Phi_approx(cS1[i, k] - S1mu[i]))/C_area_rS1;
    }
    pr[nratings] = (Phi_approx(c1[i] - S1mu[i])-Phi_approx(cS1[i, nratings-1] - S1mu[i]))/C_area_rS1;   
      
      
    // Get nI_rS2 probs (False alarm)
    pr[nratings+1] = ((1-Phi_approx(c1[i] - S1mu[i]))-(1-Phi_approx(cS2[i, 1]-S1mu[i])))/I_area_rS2;  
    for (k in 1:(nratings-2)) {
      pr[nratings+1+k] = ((1-Phi_approx(cS2[i, k] - S1mu[i]))-(1-Phi_approx(cS2[i, k+1] - S1mu[i])))/I_area_rS2;
    }
    pr[nratings*2] = (1-Phi_approx(cS2[i, nratings-1] - S1mu[i]))/I_area_rS2;  
      
      
    // Get nI_rS1 probs (Miss)
    pr[nratings*2+1] = Phi_approx(cS1[i, 1]-S2mu[i])/I_area_rS1;
    for (k in 1:(nratings-2)){
      pr[(nratings*2)+1+k] = (Phi_approx(cS1[i, k+1]-S2mu[i])-Phi_approx(cS1[i, k]-S2mu[i]))/I_area_rS1;
    }
    pr[nratings*3] = (Phi_approx(c1[i]-S2mu[i])-Phi_approx(cS1[i, nratings-1]-S2mu[i]))/I_area_rS1;
    
      
    // Get nC_rS2 probs (Hit)
    pr[(nratings*3)+1] = ((1-Phi_approx(c1[i]-S2mu[i]))-(1-Phi_approx(cS2[i, 1]-S2mu[i])))/C_area_rS2;
    for (k in 1:(nratings-2)){
      pr[(nratings*3)+1+k] = ((1-Phi_approx(cS2[i, k]-S2mu[i]))-(1-Phi_approx(cS2[i, k+1]-S2mu[i])))/C_area_rS2;
    }
    pr[nratings*4] = (1-Phi_approx(cS2[i, nratings-1]-S2mu[i]))/C_area_rS2;
      
    //Avoid underflow of probabilities
    for (j in 1:4) {
    prT[(j-1)*nratings+1:j*nratings] = normalize(0.999 * pr[(j-1)*nratings+1:j*nratings] + 0.001);
    }
      
    // Multinomial likelihood for response counts ordered as c(nR_S1,nR_S2)
    log_lik[i, 1] = multinomial_lpmf(counts[i, 1:nratings] | prT[1:nratings]); //CR
    log_lik[i, 2] = multinomial_lpmf(counts[i, (nratings+1):(nratings*2)] | prT[(nratings+1):(nratings*2)]); //False alarm
    log_lik[i, 3] = multinomial_lpmf(counts[i, ((nratings*2)+1):(nratings*3)] | prT[((nratings*2)+1):(nratings*3)]); //Miss
    log_lik[i, 4] = multinomial_lpmf(counts[i, ((nratings*3)+1):(nratings*4)] | prT[((nratings*3)+1):(nratings*4)]); //Hit
  } // end for i loop
}



model {
///// Type 2 SDT MODEL (META-d) /////

  //Type 2 priors
  meta_d ~ normal(d1,0.5);       //sampling meta-d 

  for (i in 1:N){                //sampling c2 vectors for S1 and S2
    for (j in 1:(nratings-1)) {
      cS1_p[i, j] ~ uniform(0,1); //c2 vectors (cS1_p -> cS1_raw -> cS1)
      cS2_p[i, j] ~ uniform(0,1); //c2 vectors (cS2_p -> cS2_raw -> cS2)
      }
  }
  
  target += sum(log_lik);
}

generated quantities{
  real y_pred[N, 4*nratings];
  vector[N] Mratio;
  
  for (i in 1:N){
    Mratio[i] = meta_d[i]/d1[i];
  }


  ////For posterior predictive check
  
  //Set all posterior predictions to 0 (avoids NULL values)
  
  for (i in 1:N){
    for (r in 1:(4*nratings)){
      y_pred[i,r] = -1;
    } //end of r loop
  } //end of i loop
  
  
  { //local section, this saves time and space
    for (i in 1:N){
      real C_area_rS1; // Correct rejection
      real I_area_rS1; // Miss
      real C_area_rS2; // Hit
      real I_area_rS2; // False Alarm
      vector[nratings*4] pr;
      vector[nratings*4] prT;
      
      // Calculate normalization constants
      C_area_rS1 = Phi_approx(c1[i] - S1mu[i]);
      I_area_rS1 = Phi_approx(c1[i] - S2mu[i]);
      C_area_rS2 = 1-Phi_approx(c1[i] - S2mu[i]);
      I_area_rS2 = 1-Phi_approx(c1[i] - S1mu[i]);
    
      // Get nC_rS1 probs (Correct Rejection)
      pr[1] = Phi_approx(cS1[i, 1] - S1mu[i])/C_area_rS1;  
      for (k in 1:(nratings-2)) {
        pr[k+1] = (Phi_approx(cS1[i, k+1] - S1mu[i])-Phi_approx(cS1[i, k] - S1mu[i]))/C_area_rS1;
      }
      pr[nratings] = (Phi_approx(c1[i] - S1mu[i])-Phi_approx(cS1[i, nratings-1] - S1mu[i]))/C_area_rS1;   
      // Get nI_rS2 probs (False alarm)
      pr[nratings+1] = ((1-Phi_approx(c1[i] - S1mu[i]))-(1-Phi_approx(cS2[i, 1]-S1mu[i])))/I_area_rS2;  
      for (k in 1:(nratings-2)) {
        pr[nratings+1+k] = ((1-Phi_approx(cS2[i, k] - S1mu[i]))-(1-Phi_approx(cS2[i, k+1] - S1mu[i])))/I_area_rS2;
      }
      pr[nratings*2] = (1-Phi_approx(cS2[i, nratings-1] - S1mu[i]))/I_area_rS2;  
      // Get nI_rS1 probs (Miss)
      pr[nratings*2+1] = Phi_approx(cS1[i, 1]-S2mu[i])/I_area_rS1;
      for (k in 1:(nratings-2)){
        pr[(nratings*2)+1+k] = (Phi_approx(cS1[i, k+1]-S2mu[i])-Phi_approx(cS1[i, k]-S2mu[i]))/I_area_rS1;
      }
      pr[nratings*3] = (Phi_approx(c1[i]-S2mu[i])-Phi_approx(cS1[i, nratings-1]-S2mu[i]))/I_area_rS1;
      // Get nC_rS2 probs (Hit)
      pr[(nratings*3)+1] = ((1-Phi_approx(c1[i]-S2mu[i]))-(1-Phi_approx(cS2[i, 1]-S2mu[i])))/C_area_rS2;
      for (k in 1:(nratings-2)){
        pr[(nratings*3)+1+k] = ((1-Phi_approx(cS2[i, k]-S2mu[i]))-(1-Phi_approx(cS2[i, k+1]-S2mu[i])))/C_area_rS2;
      }
      pr[nratings*4] = (1-Phi_approx(cS2[i, nratings-1]-S2mu[i]))/C_area_rS2;
      //Avoid underflow of probabilities
      //Avoid underflow of probabilities
      for (j in 1:4) {
        prT[(j-1)*nratings+1:j*nratings] = normalize(0.999 * pr[(j-1)*nratings+1:j*nratings] + 0.001);
      } //end of j loop

      y_pred[i, 1:nratings] = multinomial_rng(prT[1:nratings], CR[i]); //CR
      y_pred[i, (nratings+1):(nratings*2)] = multinomial_rng(prT[(nratings+1):(nratings*2)], FA[i]); //False alarm
      y_pred[i, (nratings*2+1):(nratings*3)] = multinomial_rng(prT[(nratings*2+1):(nratings*3)], M[i]); // Miss
      y_pred[i, (nratings*3+1):(nratings*4)] = multinomial_rng(prT[(nratings*3+1):(nratings*4)], H[i]); // Hit rate

    } //end of i loop
  } //end of local section
}
