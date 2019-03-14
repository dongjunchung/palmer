#if !defined(_Main_H)
#define _Main_H

#include <RcppArmadillo.h>
#include "3_Param.h"

class CMain {
	
  public:
		
    CMain(arma::mat Y_);
	  ~CMain(); //destructor

    void Initialize() ; void check_random_generate() ; 
    void Iterate(); void Run(int iter);
    void SetMsgLevel(int level); int GetMsgLevel();
    void clearE_ijt() ; 
    // void SetOptionEit(bool OptionEit_); bool GetOptionEit();
    
		// Initial values // 
    void SetE_mat(arma::mat E_mat_); arma::mat GetE_mat();
		void Setmu_vec(arma::vec mu_vec_); arma::vec Getmu_vec();
		void Setsig2_vec(arma::vec sig2_vec_); arma::vec Getsig2_vec();
    void SetBeta(arma::mat Beta_); arma::mat GetBeta();
    void SetG_mat(arma::mat G_mat); arma::mat GetG_mat();
				
		// Hyperparameters	
		void Settheta_mu(double theta_mu_) ; double Gettheta_mu() ;  
		void Settau2_mu(double tau2_mu_) ; double Gettau2_mu() ;
		void Seta_sigma(double a_sigma_) ; double Geta_sigma() ; 
		void Setb_sigma(double b_sigma_) ; double Getb_sigma() ; 
		void Settheta_alpha(double theta_alpha_) ; double Gettheta_alpha(); 
		void Settau2_alpha(double tau2_alpha_) ; double Gettau2_alpha() ;
		void Setstepsize_alpha(double stepsize_alpha_) ; double Getstepsize_alpha() ;
		void Seta_beta(double a_beta_) ; double Geta_beta() ;
		void Setb_beta(double b_beta_) ; double Getb_beta() ;
		void Setstepsize_beta(double stepsize_beta_) ; double Getstepsize_beta() ;
		void Seta_betaG(double a_betaG_) ; double Geta_betaG() ;
		void Setb_betaG(double b_betaG_) ; double Getb_betaG() ;
		void Setthreshold_on(double threshold_on_) ; double Getthreshold_on() ;
		
		// Print data or result
		arma::mat GetY();
		arma::vec GetAccept() ; arma::vec GetAccProb();
		double GetnormC() ;
		double GetlogPost() ;
		double Getloglikelihood() ;
		arma::cube Getsum_E_ijt() ;
		
  private:
    
    CData Data;
    CParam Param;
    int IterCount;

};

#endif  //_CMain_H
