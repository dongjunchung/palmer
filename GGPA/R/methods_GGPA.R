
# generic methods for "GGPA" class

# GGPA model fit summary

setMethod(
    f="show",
    signature="GGPA",
    definition=function( object ) {    
    
    # summary of GGPA fit
		
    # constants
		
		nBin <- nrow(object@gwasPval)
		nGWAS <- ncol(object@gwasPval)
		
		# estimates
		
		est_mu_vec = sd_mu_vec = rep(0,nGWAS)
    est_sigma1 = sd_sigma1 = rep(0,nGWAS)
		
    for (i in 1:nGWAS){
      est_mu_vec[i] = mean(object@fit$mu[,i])
      sd_mu_vec[i] = sd(object@fit$mu[,i])
      est_sigma1[i] = mean(object@fit$sigma[,i])
      sd_sigma1[i] = sd(object@fit$sigma[,i])
    }
    
    est_mean_E = apply(object@fit$mean_E,2,mean)
    sd_mean_E = apply(object@fit$mean_E,2,sd)
		
    MU = round(cbind(est_mu_vec,sd_mu_vec),2)
    rownames(MU) = colnames(object@gwasPval)
    colnames(MU) <- c( "estimate", "SE" )
    
    SIGMA = round(cbind(est_sigma1,sd_sigma1),2)
    rownames(SIGMA) = colnames(object@gwasPval) 
    colnames(SIGMA) <- c( "estimate", "SE" )
    
    EMAT = round(cbind(est_mean_E,sd_mean_E),2)
    rownames(EMAT) = colnames(object@gwasPval) 
    colnames(EMAT) <- c( "estimate", "SE" )
		
		# output
        
    cat( "Summary: GGPA model fitting results (class: GGPA)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Data summary:\n" )
    cat( "\tNumber of GWAS data: ", nGWAS , "\n", sep="" )
		cat( "\tNumber of SNPs: ", nBin , "\n", sep="" )
		cat( "mu\n", sep="" )
		print(MU)
		cat( "sigma\n", sep="" )
		print(SIGMA)
		cat( "Proportion of associated SNPs\n", sep="" )
		print(EMAT)
    cat( "--------------------------------------------------\n" )
  }
)

# phenotype graph

setMethod(
  f="plot",
  signature=c("GGPA","missing"),
  definition=function( x, y, pCutoff = 0.5, betaCI = 0.95, ... ) {   
     
    P_hat_ij <- x@summary$P_hat_ij
    draw_beta <- x@fit$beta
     
    # calculate posterior probabilities
    
    Names = dimnames(P_hat_ij)[[1]]
    n_pheno = dim(draw_beta)[[2]]
    
    P_hat = P_lb_beta = matrix( 0, n_pheno, n_pheno )
    for (i in 1:n_pheno){
    	for (j in 1:n_pheno){
    		P_hat[i,j] = mean( draw_beta[,i,j] > 0 )
    		P_lb_beta[i,j] = quantile( draw_beta[,i,j], probs = ( 1 - betaCI ) / 2 )
    	}
    }
    dimnames(P_hat)[[1]] = dimnames(P_lb_beta)[[1]] = Names
    dimnames(P_hat)[[2]] = dimnames(P_lb_beta)[[2]] = Names
    
    # graph construction	
    
    adjmat <- round(P_hat,2) > pCutoff & P_lb_beta > 0
    rownames(adjmat) <- colnames(adjmat) <- Names
    
    am.graph <- new( "graphAM", adjMat=adjmat, edgemode="undirected" )
    am.graph
    plot( am.graph )
  }
)

# local FDR

setMethod(
    f="fdr",
    signature="GGPA",
    definition=function( object, i=NULL, j=NULL ) {
        # return marginal FDR
		
    # constants
		
		nBin <- nrow(object@gwasPval)
		nGWAS <- ncol(object@gwasPval)
		
		if ( is.null(i) & is.null(j) ) {
		  #message( "Info: Marginal local FDR matrix is returned." )
		  
  		fdrmat <- matrix( NA, nBin, nGWAS )
  		colnames(fdrmat) <- colnames(object@gwasPval)
  		
  		for (i_phen in 1:nGWAS){
    		Prob_e_ijt = object@summary$Sum_E_ijt[i_phen,i_phen,] / length(object@fit$loglik)
        fdrmat[,i_phen] <- 1 - Prob_e_ijt
  		}
		} else if ( !is.null(i) & !is.null(j) ) {
		  #message( "Info: Local FDR vector for specified i & j pair is returned." )
		  
		  fdrmat <- rep( NA, nBin )
		  
		  Prob_e_ijt = object@summary$Sum_E_ijt[i,j,] / length(object@fit$loglik)
      fdrmat <- 1 - Prob_e_ijt
		} else {
		  stop( "Both of i and j should be either NULL or numeric!" )
		}
    
    return(fdrmat)
  }
)

# parameter estimates

setMethod(
    f="estimates",
    signature="GGPA",
    definition=function( object, ... ) {
        # return parameter estimates
		
		return(object@summary)
  }
)
