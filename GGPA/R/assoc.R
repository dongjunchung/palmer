
setMethod(
    f="assoc",
    signature="GGPA",
    definition=function( object, FDR=0.05, fdrControl="global", i=NULL, j=NULL ) {
		
		# check arguments
		
		if ( fdrControl != "global" & fdrControl != "local" ) {
			stop( "Invalid value for 'fdrControl' argument! It should be either 'global' or 'local'." )
		}
		
		if ( FDR < 0 | FDR > 1 ) {
			stop( "Invalid value for 'FDR' argument! It should be between zero and one." )
		}
		
		# association mapping
		
		fdrmat <- fdr( object, i, j )
    
    if ( is.null(i) & is.null(j) ) {
  		# based on marginal FDR
      
      message( "Info: Association mapping for each phenotype." )
      
      amat <- matrix( 0, nrow(fdrmat), ncol(fdrmat) )
  		
  		if ( fdrControl == "local" ) {
  			# local FDR control
  			
  			message( "Info: Association mapping based on the local FDR control at level ", FDR, "." )
  			
  			amat[ fdrmat <= FDR ] <- 1
  		} else if ( fdrControl == "global" ) {
  			# global FDR control
  			
  			message( "Info: Association mapping based on the global FDR control at level ", FDR, "." )
  		
  			# direct approach for FDR control
  			
  			for ( j_pheno in 1:ncol(amat) ) {
  				pp <- fdrmat[,j_pheno]
  				pp.ordered <- sort(pp)
  				pp.cum <- cumsum( pp.ordered ) / c(1:length(pp))
  				cutoff <- max( pp.ordered[ pp.cum <= FDR ] )
  				amat[ pp <= cutoff, j_pheno ] <- 1
  			}  			
  		}
    } else if ( !is.null(i) & !is.null(j) ) {
      # based on FDR of interest
      
      message( "Info: Association mapping for specified i & j phenotype pair." )
      
      amat <- rep( 0, length(fdrmat) )
      
      if ( fdrControl == "local" ) {
        # local FDR control
        
        message( "Info: Association mapping based on the local FDR control at level ", FDR, "." )
        
        amat[ fdrmat <= FDR ] <- 1
      } else if ( fdrControl == "global" ) {
        # global FDR control
        
        message( "Info: Association mapping based on the global FDR control at level ", FDR, "." )
        
        # direct approach for FDR control
        
        pp <- fdrmat
        pp.ordered <- sort(pp)
        pp.cum <- cumsum( pp.ordered ) / c(1:length(pp))
        cutoff <- max( pp.ordered[ pp.cum <= FDR ] )
        amat[ pp <= cutoff ] <- 1
        
      }      
    }
		
		return(amat)
	}
)
