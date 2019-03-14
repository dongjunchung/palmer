#include "3_Param.h" 

#define LOG_2_PI 1.83787706640935

CData::CData(){ Debug = false ; }
CData::~CData(){ } //Destructor

CParam::CParam(){
  accept_prob_vec = arma::zeros<arma::vec>(9) ; 
  is_accept_vec = arma::zeros<arma::vec>(9) ; 
}
CParam::~CParam(){ }

void CParam::Initialize(CData &Data){
  n_pheno = Data.Y.n_rows ; n_SNP = Data.Y.n_cols ; // constants but stored in CParam for convienence
  Data.logY = arma::zeros<arma::mat>(n_pheno,n_SNP) ; 
		// Note 1. This is used to calculate sum of log y_it when e_it=1. 
		//      2. (Weak signal) Because e_it=0 with y_it <= 0 by definition (model assumption),
		//          we will not use log y_it for e_it=0 and store zero here. 
		//			3. (Strong signal) Also, we put e_it=1 when y_it > threshold    // V 1.3.1 
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++ ){
    for (int i_SNP=0; i_SNP<n_SNP; i_SNP++ ){
      if ( Data.Y(i_pheno,i_SNP) > 0){ // V 2.0.2
        Data.logY(i_pheno,i_SNP) = log(Data.Y(i_pheno,i_SNP)) ;
				if ( Data.Y(i_pheno,i_SNP) > Data.threshold_on ) E_mat(i_pheno,i_SNP) = 1 ; 
      } else {
      	E_mat(i_pheno,i_SNP) = 0 ;
      }
    }
  } 
  
  normC = normC_fn(Beta, Data) ;	// Ver_1_4_1
  if ( normC < 0 ){
    Rcpp::stop("The initialized normC has a negative value.") ;
  }
  Data.msg_level = 0; //0 errors only; 1: errors and warnings; 2: errors, warnings and information
  
  sum_E_ijt = arma::zeros<arma::cube>(n_pheno,n_pheno,n_SNP) ; 
  
  is_initialized = 1 ; 
}

void CParam::check_random_generate(CData &Data){
  RandVec = Rcpp::rnorm(2,2,1000) ;
  std::cout << RandVec(0) << "  " << RandVec(1) << std::endl ; 
  // std::cout << pow(10,-3) << std::endl ; 
}

void CParam::iterate(int iter, CData &Data, int n_simul) {
  if ( is_initialized == 1 ){
    logPost = 0.0 ; loglikelihood = 0.0 ; 
    S1_e_it(Data) ;
    S2_mu_i(Data) ;  
    S3_sig2_i(Data) ;
    S4_alpha_i(Data) ; 
  	S5_beta_ij(Data) ; 
  	S6_G_beta_ij(Data) ;  
  	store_Eit(Data) ; 
  	// if ( Data.OptionEit==TRUE ) store_Eit(Data) ; 
  } else {
    Rcpp::stop("Need To Run model$Initialize()") ; 
  }
} 

void CParam::clearE_ijt( ){
  sum_E_ijt = arma::zeros<arma::cube>(n_pheno,n_pheno,n_SNP) ; 
}

///////////
void CParam::S1_e_it(CData &Data) {
	
  for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
    for (int i=0; i<n_pheno; i++){
      if ( ( Data.Y(i,i_SNP) > 0 ) & ( Data.Y(i,i_SNP) <= Data.threshold_on ) ){ 
				// V 1.3.1 -> other case: e_it is fixed in CParam::Initialize, 
				//             i.e., e_it=0 if y_it<=0 and 1 if y_it>=Data.threshold_on
        double unnorm_logprob0 = 0.0 ;
        double unnorm_logprob1 = Beta(i,i) ;
        for (int j=0; j<n_pheno; j++){
          if (G_mat(i,j)==1) unnorm_logprob1 = unnorm_logprob1 + Beta(i,j) * E_mat(j,i_SNP) ; 
        }
				  // Note: Not count G_mat(i,j)=0 or G_mat(i,j)=9,i.e., diagonal
        unnorm_logprob0 = unnorm_logprob0 + R::dnorm(Data.Y(i,i_SNP),0,1,1) ; // log = T
        unnorm_logprob1 = unnorm_logprob1 + R::dlnorm(Data.Y(i,i_SNP),mu_vec(i),sqrt(sig2_vec(i)),1) ; // log = T
        double prob_e_it = 1.0 / ( 1.0 + exp(unnorm_logprob0-unnorm_logprob1) ) ; 
        	// Note: p1 = f1/(f1+f0) = 1 / (1+f0/f1) = 1 / (1+exp(log f0 - log f1)) 
				RandVec = Rcpp::runif(1,0,1) ;
        if ( prob_e_it >= RandVec(0) ){
          E_mat(i,i_SNP) = 1 ; 
        } else {
          E_mat(i,i_SNP) = 0 ; 
        }
      }         
      // To check loglikelihood
      if (E_mat(i,i_SNP)==1){
        loglikelihood = loglikelihood + R::dlnorm(Data.Y(i,i_SNP),mu_vec(i),sqrt(sig2_vec(i)),1)  ;
      } else {
        loglikelihood = loglikelihood + R::dnorm(Data.Y(i,i_SNP),0,1,1) ;
      }
		} 
  }
	
  is_accept_vec(0) = 1 ; 
} 

//////////
void CParam::S2_mu_i(CData &Data) {

	for (int i=0; i<n_pheno; i++){
		double sig2_i = sig2_vec(i) ;
		arma::vec e_i = E_mat.row(i).t() ; 
		double n_i = sum(e_i) ; 
		arma::vec logy_i = Data.logY.row(i).t() ;  
		arma::vec e1_logy = logy_i % e_i ; 
		double sum_e1_logy = sum(e1_logy) ; 
			// Note: % Schur product: element-wise multiplication of two objects
			//         sum of log_y_it with e_it=1
		double mean_star = (sig2_i * Data.theta_mu + Data.tau2_mu * sum_e1_logy) / (sig2_i + Data.tau2_mu * n_i) ; 
		double var_star = (sig2_i * Data.tau2_mu)/(sig2_i + Data.tau2_mu * n_i) ; 
		RandVec = Rcpp::rnorm(1, mean_star, sqrt(var_star)) ; 
		mu_vec(i) = RandVec(0) ;
  } 
	
	is_accept_vec(1) = 1 ; 
} 

///////////
void CParam::S3_sig2_i(CData &Data) {

	for (int i=0; i<n_pheno; i++){
		arma::vec mu_i_onevec(n_SNP) ; mu_i_onevec.fill(mu_vec(i)) ; 
		arma::vec e_i = E_mat.row(i).t() ; 
		double n_i = sum(e_i) ; 
		arma::vec logy_i = Data.logY.row(i).t() ;  
		arma::vec logy_minus_mu = logy_i - mu_i_onevec ;  
		arma::vec e1_logy_minus_mu = logy_minus_mu % e_i ;
		arma::vec sum_e1_squares = e1_logy_minus_mu.t() * e1_logy_minus_mu ; 
		double a_star = Data.a_sigma + 0.5 * n_i ; 
		double b_star = Data.b_sigma + 0.5 * sum_e1_squares(0) ; 
		sig2_vec(i) = rinvgamma(a_star, b_star) ; 
  } 

	is_accept_vec(2) = 1 ;
}

///////////
void CParam::S4_alpha_i(CData &Data) {
  
  is_accept_vec(3) = 0 ; 
  // if ( normC < 0 ) normC = normC_fn(Beta, Data) ;	// Ver_1_4_1
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){  // V_1_3_5 ; Update each alpha_i
    arma::mat Beta_prop = Beta ;
    RandVec = Rcpp::rnorm(1, Beta(i_pheno,i_pheno), Data.stepsize_alpha ) ;
    Beta_prop(i_pheno,i_pheno) = RandVec(0) ; 
    double logP_numer = R::dnorm(Beta_prop(i_pheno,i_pheno), Data.theta_alpha, sqrt(Data.tau2_alpha),1) ; // log=TRUE
    double logP_denom = R::dnorm(Beta(i_pheno,i_pheno), Data.theta_alpha, sqrt(Data.tau2_alpha),1) ; // log=TRUE
    double normC_prop = normC_fn(Beta_prop, Data) ;
    for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
      double e_it = E_mat(i_pheno,i_SNP) ;
      logP_numer = logP_numer + log( normC ) + Beta_prop(i_pheno,i_pheno) * e_it  ; 
      logP_denom = logP_denom + log( normC_prop ) + Beta(i_pheno,i_pheno) * e_it ;
    }
    double accept_prob = exp( logP_numer - logP_denom ) ; 
    accept_prob_vec(3) = accept_prob ;  
    RandVec = Rcpp::runif(1, 0, 1) ; 
    if ( accept_prob >= RandVec(0) ){
      Beta = Beta_prop ; normC = normC_prop ; 
      is_accept_vec(3) = is_accept_vec(3) + 1 ; 
    } 
  } // for (int i=0; i<n_pheno; i++){

  is_accept_vec(3) = 1.0 / n_pheno * is_accept_vec(3) ;
} 

///////////
void CParam::S5_beta_ij(CData &Data) {
  
	is_accept_vec(4) = 0 ; 
	int w_cur = 0 ; // just for calculation of is_accept_vec( ) 
	for (int i=0; i<(n_pheno-1); i++){
		for (int j=(i+1); j<n_pheno; j++){
			if ( G_mat(i,j)==1 ){
				w_cur ++ ; 
			  arma::mat Beta_prop = Beta ;
			  
				double temp_beta_prop = 0.0 ; 
				while ( temp_beta_prop <= 0.0 ){
					temp_beta_prop = rtruncNorm_uppertail_fn(Beta(i,j), Data.stepsize_beta, 0) ;
				}
				Beta_prop(i,j) = temp_beta_prop ; 
				Beta_prop(j,i) = temp_beta_prop ;
				double logP_numer = R::dgamma(Beta_prop(i,j),Data.a_beta,(1.0/Data.b_beta),1) ;   // (shape,scale,log), i.e., b_beta = rate // V_1_3_5
			  double logP_denom = R::dgamma(Beta(i,j),Data.a_beta,(1.0/Data.b_beta),1) ;         // 
        double normC_prop = normC_fn(Beta_prop, Data) ;
			  for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
			    double e_it = E_mat(i,i_SNP) ;
			    double e_jt = E_mat(j,i_SNP) ; 
			    logP_numer = logP_numer + log(normC) + Beta_prop(i,j) * e_it * e_jt ; 
			    logP_denom = logP_denom + log(normC_prop) + Beta(i,j) * e_it * e_jt ; 
			  }
			  double logQ_numer = log( dtruncnorm_uppertail_fn(Beta(i,j), 0, Beta_prop(i,j), Data.stepsize_beta) ) ; // CHECK
			  double logQ_denom = log( dtruncnorm_uppertail_fn(Beta_prop(i,j), 0, Beta(i,j), Data.stepsize_beta) ) ;
				double accept_prob = exp( logP_numer - logP_denom + logQ_numer - logQ_denom ) ; 
			  accept_prob_vec(4) = accept_prob ;  
        RandVec = Rcpp::runif(1, 0, 1) ; 
			  if ( accept_prob >= RandVec(0) ){
			    Beta = Beta_prop ; normC = normC_prop ; 
					is_accept_vec(4) = is_accept_vec(4) + 1 ; 
			  } 
			}
		}
	}
	
	if ( w_cur == 0 ){
	  is_accept_vec(4) = 0.2 ; 
	} else {
	  is_accept_vec(4) = 1.0 / w_cur * is_accept_vec(4) ; 
	}
} 

///////////
void CParam::S6_G_beta_ij(CData &Data) {
  
  arma::mat G_prop = G_mat ; arma::mat Beta_prop = Beta ; 
  double logQ_numer, logQ_denom ; 
  double logP_numer, logP_denom ;
  double normC_prop ;
  
  // Step 1
  int w_cur = 0 ; int w_max = 0 ; int w_prop ; 
  for (int i=0; i<(n_pheno-1); i++){
    for (int j=(i+1); j<n_pheno; j++){
      w_max ++ ; w_cur = w_cur + G_mat(i,j) ; 
    }
  }
  if ( w_cur==0 ){ w_prop = 1 ; logQ_numer = log(0.5) ; logQ_denom = log(1.0) ; } // q(w|w^q) // q(w^q|w)
  if ( w_cur==w_max ){ w_prop = w_max - 1 ; logQ_numer = log(0.5) ;  logQ_denom = log(1.0) ; } // q(w|w^q) // q(w^q|w)
  if ( (w_cur > 0) && (w_cur < w_max) ){
    RandVec = Rcpp::runif(1, 0, 1) ; 
    if ( RandVec(0) < 0.5 ){ w_prop = w_cur + 1 ; } else { w_prop = w_cur - 1 ; }
    logQ_denom = log(0.5) ; // q(w^q|w)
    if ( (w_prop==0) || (w_prop==w_max) ){ logQ_numer = log(1.0) ; } else { logQ_numer = log(0.5) ; } // q(w|w^q)
  }
  
  // Step 2
  if ( w_prop > w_cur  ){ // step 2-a
    int id_added = rDiscrete(w_max-w_cur) + 1 ; // 1 ~ total. of empty edges
    int count_empty_edge = 0 ;
    for (int i=0; i<(n_pheno-1); i++){
      for (int j=(i+1); j<n_pheno; j++){
        if ( G_mat(i,j)==0 ){
          count_empty_edge++;
          if ( id_added==count_empty_edge ){
            G_prop(i,j) = 1 ; G_prop(j,i) = G_prop(i,j) ; 
            RandVec = Rcpp::rgamma(1, Data.a_betaG, 1.0/Data.b_betaG) ; // q(beta_ij^q)
            Beta_prop(i,j) = RandVec(0) ; 
            Beta_prop(j,i) = Beta_prop(i,j) ; 
            id_added = -9 ; 
            
            logP_numer = R::dgamma(Beta_prop(i,j),Data.a_beta,(1.0/Data.b_beta),1) ; // f(beta_ij^q|G_ij^q)
            logP_denom = 0 ; // f(beta_ij|G_ij) cancelled with q(beta_ij)
            normC_prop = normC_fn(Beta_prop, Data) ;
            for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
              double e_it = E_mat(i,i_SNP) ;
              double e_jt = E_mat(j,i_SNP) ; 
              logP_numer = logP_numer - log(normC_prop) + Beta_prop(i,j) * e_it * e_jt ; // f(e_t|alpha,beta^q,G^q)
              logP_denom = logP_denom - log(normC) + Beta(i,j) * e_it * e_jt ; // f(e_t|alpha,beta,G)
            }
            logQ_numer = logQ_numer + 0 ; // q(beta_ij) cancelled with f(beta_ij|G_ij)  
            logQ_denom = logQ_denom + R::dgamma(Beta_prop(i,j),Data.a_betaG,(1.0/Data.b_betaG),1) ; // q(beta_ij^q)
          }
        } 
      } 
    }
    logQ_numer = logQ_numer - log(w_cur) ; // q(G|G^q,w)
    logQ_denom = logQ_denom - log(w_max-w_cur)  ; // q(G^q|G,w^q)
  } else { // step 2-b
    int id_deleted = rDiscrete(w_cur)+1 ; // 1 ~ total no. of connected edges
    int count_connected_edge = 0 ;
    for (int i=0; i<(n_pheno-1); i++){
      for (int j=(i+1); j<n_pheno; j++){
        if ( G_mat(i,j)==1 ){
          count_connected_edge++;
          if ( id_deleted==count_connected_edge ){
            G_prop(i,j) = 0 ; G_prop(j,i) = G_prop(i,j) ; 
            Beta_prop(i,j) = 0 ; // q(beta_ij^q)
            Beta_prop(j,i) = Beta_prop(i,j) ; 
            id_deleted = -9 ; 
            
            logP_numer = 0 ; // f(beta_ij^q|G_ij^q) cancelled with q(beta_ij^q)
            logP_denom = R::dgamma(Beta(i,j),Data.a_beta,(1.0/Data.b_beta),1) ; // f(beta_ij|G_ij)
            normC_prop = normC_fn(Beta_prop, Data) ;
            for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
              double e_it = E_mat(i,i_SNP) ;
              double e_jt = E_mat(j,i_SNP) ; 
              logP_numer = logP_numer - log(normC_prop) + Beta_prop(i,j) * e_it * e_jt ; // f(e_t|alpha,beta^q,G^q)
              logP_denom = logP_denom - log(normC) + Beta(i,j) * e_it * e_jt ; // f(e_t|alpha,beta,G)
            }
            logQ_numer = logQ_numer + R::dgamma(Beta(i,j),Data.a_betaG,(1.0/Data.b_betaG),1) ; // q(beta_ij) 
            logQ_denom = logQ_denom + 0 ; // q(beta_ij^q) cancelled with f(beta_ij^q|G_ij^q) 
          }
        } 
      }
    }
    logQ_numer = logQ_numer - log(w_max-w_cur) ; // q(G|G^q,w)
    logQ_denom = logQ_denom - log(w_cur) ; // q(G^q|G,w^q)
  }
  if ( w_prop > 0 ) logP_numer = logP_numer - log(w_prop) ; // f(G^q) // if w_prop=0, let 1/w_prop = 1, so that log(w_prop) = 0
  if ( w_cur > 0 ) logP_denom = logP_denom - log(w_cur) ; // f(G) // if w_cur=0, let 1/w_cur = 1, so that log(w_cur) = 0

  // Step 3  
  double accept_prob = exp( logP_numer - logP_denom + logQ_numer - logQ_denom ) ; 
  accept_prob_vec(5) = accept_prob ;  
  RandVec = Rcpp::runif(1, 0, 1) ; 
  if ( accept_prob >= RandVec(0) ){
    Beta = Beta_prop ; G_mat = G_prop ; normC = normC_prop ; 
    is_accept_vec(5) = 1 ; 
  } else {
    is_accept_vec(5) = 0 ; 
  }
} 

// ///////////
void CParam::store_Eit(CData &Data) {
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
    for (int j_pheno=0; j_pheno<n_pheno; j_pheno++){
      for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
        sum_E_ijt(i_pheno,j_pheno,i_SNP) = sum_E_ijt(i_pheno,j_pheno,i_SNP) + E_mat(i_pheno,i_SNP) * E_mat(j_pheno,i_SNP) ; 
      }   
    }   
  }
}

// void CParam::store_Eit(CData &Data) {
//   
//   FILE *FILE_E_it; 
//   
//   std::string Name1 = "E_it_1" ; 
//   Name1.append(".dat") ; 
//   
//   FILE_E_it = fopen(Name1.c_str(),"a");
//   for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
//     fprintf(FILE_E_it, "%.0f ",E_mat(0,i_SNP));
//   }
//   fprintf(FILE_E_it, "\n");
//   fclose(FILE_E_it);
//   
//   std::string Name2 = "E_it_2" ; 
//   Name2.append(".dat") ; 
//   
//   FILE_E_it = fopen(Name2.c_str(),"a");
//   for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
//     fprintf(FILE_E_it, "%.0f ",E_mat(0,i_SNP));
//   }
//   fprintf(FILE_E_it, "\n");
//   fclose(FILE_E_it);
//   
// }
  
//////////////////////////////////////
// Function

double CParam::normC_fn(arma::mat Beta_input, CData &Data){
  
    double normC_temp = 1.0 ; // exp(0) 
    arma::vec e_temp(n_pheno) ; 
    int n_pow_pheno = pow(2,n_pheno) ;
    
    for (int i_pow=1; i_pow<n_pow_pheno; i_pow++){ // exclude (0 0 0)
      
      int temp_no = i_pow ; 
      for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
        int temp_denominator = pow(2,(n_pheno-i_pheno-1)) ; 
        if ( floor(temp_no/temp_denominator) == 1 ){
          e_temp(i_pheno) = 1 ; 
          temp_no = temp_no - temp_denominator ; 
        } else {
          e_temp(i_pheno) = 0 ; 
        }
      } // for i_pheno
      
      double temp_beta_sum = 0.0 ; 
      for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
        if (e_temp(i_pheno)==1){
          temp_beta_sum = temp_beta_sum + Beta_input(i_pheno,i_pheno)  ; 
          if (i_pheno<(n_pheno-1)){
            for (int j_pheno=(i_pheno+1); j_pheno<n_pheno; j_pheno++){
              if (e_temp(j_pheno)==1){
                temp_beta_sum = temp_beta_sum + Beta_input(i_pheno,j_pheno) ; 
              } // if (e_temp(j_pheno)==1)
            } // for j_pheno 
          } // if (i_pheno<(n_pheno-1))
        } // if (e_temp(i_pheno)==1)
      } // for i_pheno
      
      normC_temp = normC_temp + exp(temp_beta_sum) ; 
      
    } // for i_pow
    
    return normC_temp ;
} 

//////////////////////////////////////
// Distribution

double CParam::rinvgamma(double alpha, double beta){
  RandVec = Rcpp::rgamma(1, alpha, (1.0/beta) ) ; 
  return ( 1.0 / RandVec(0) ) ; 
  // Note that b in Rcpp::rgamma(a,b) is scale, 
  // i.e., its mean is ab, NOT a/b  
  // If X ~ Gamma(a,b), then 1/X ~ IG(a,1/b)
}

double CParam::rtruncNorm_lowertail_fn(double mean, double sd, double upperlimit){
  RandVec = Rcpp::runif(1,0,1) ; 
  double std_upperlimit = (upperlimit-mean) / sd ; 
  double temp_pnorm = R::pnorm(std_upperlimit,0,1,1,0) * RandVec(0) ; // x,mu,sigma,lt,lg
  double std_x = R::qnorm(temp_pnorm,0,1,1,0) ;  // p,mu,sigma,lt,lg
  return( mean + std_x * sd ) ;
}

double CParam::rtruncNorm_uppertail_fn(double mean, double sd, double lowerlimit){
  RandVec = Rcpp::runif(1,0,1) ; 
  double std_lowerlimit = (lowerlimit-mean) / sd ; 
  double temp_pnorm = (1.0-R::pnorm(std_lowerlimit,0,1,1,0)) * RandVec(0) ; // x,mu,sigma,lt,lg
  double std_x = -1.0 * R::qnorm(temp_pnorm,0,1,1,0) ;  // p,mu,sigma,lt,lg
  return( mean + std_x * sd ) ;
}

double CParam::dtruncnorm_uppertail_fn(double x, double lowerlimit, double mu, double sig){
	double alpha = (lowerlimit-mu) / sig ; 
	double z = (x-mu) / sig ;
	double res = R::dnorm(z,0,1,0) / (sig * ( 1.0-R::pnorm(alpha,0,1,1,0) ) ) ; 	
  return ( res ) ;
	
}

int CParam::rDiscrete(int max_no){
  // generate an integer from 0 to (max_no-1)
  double max_no_ = 1.0 * max_no ; 
  RandVec = Rcpp::runif(1, 0.0, max_no_) ; 
  return ( floor(RandVec(0)) ) ; 
}
