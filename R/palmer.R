#' PALMER
#'
#' A constrained biclustering Algoritm to improve pathway annotation bssed on the biomedical literature mining
#'
#' @export
#'
#' @param X matrix of binary values with n sample by p variable
#' @param path list. Gene list for the prior information (default: NULL)
#' @param K The number of gene clusters (default: 3)
#' @param L The number of GO term clusters (default: 3)
#' @param B Number of bootstrapping (default: 1000)
#'
#' @examples
#' data(sdata)
#' data(pathway)
#'
#' fit.palmer <- palmer(X=sdata,path=pathway,K=2,L=3,B=100)
#'
#' fit.palmer
#' predict(fit.palmer)
#' plot(fit.palmer)


palmer <- function(X,path=NULL,K=3,L=3,B=1000) {

  em <- function(x,theta,pi0,n,pathway=NULL) {

	zik.new <- function(x,theta,pi0,n,path=pathway){
	  I <- dim(x)[1]
      J <- dim(x)[2]
      L <- dim(theta)[2]
      K <- dim(theta)[1]
      pathn <- length(unique(path))

      pz <- array(,c(I,K))
      for (i in 1:I){
        pz1 <- NULL

        for (k in 1:K){
          temp <- NULL
          for (l in 1:L){
            temp[l] <- (theta[k,l]^x[i,l])*(1-theta[k,l])^(n[l]-x[i,l])
          }
          pz1[k] <- pi0[k]*prod(temp)
        }

        if (all(pz1==0)) {
          pz[i,] <- pz1} else {
            pz[i,] <- pz1/sum(pz1)
          }

        if (!is.null(path)){
          for (k in 1:(pathn-1)){
            if (path[i]==k) {
              if (pz[i,k]<1) {
                pz[i,k] <- 1
                pz[i,-k] <- (1-pz[i,k])/(pathn-1)
              }}}}
      }
      pz
    }


    pi.new <- function(zik) apply(zik,2,sum)/sum(zik)

    theta.new <- function(x,zik,n){
      I <- dim(x)[1]
      J <- dim(x)[2]
      L <- length(n)
      K <- dim(zik)[2]

      the <- array(,c(K,L))
      for (k in 1:K){
        for (l in 1:L){
          temp <- 0
          for (i in 1:I){
            temp[i] <- (zik[i,k]*x[i,l])/(sum(zik[,k])*n[l])
          }
          the[k,l] <- sum(temp)
        }
      }
      the[the<0.01] <- 0.01; the[the>0.99] <- 0.99
      the
    }

    llike <- function(x,theta,pi0,n){
      I <- dim(x)[1]
      J <- dim(x)[2]
      L <- dim(theta)[2]
      K <- dim(theta)[1]

      like <- NULL
      for (i in 1:I){
        pz1 <- NULL
        for (k in 1:K){
          temp <- NULL
          for (l in 1:L){
            temp[l] <- (theta[k,l]^x[i,l])*(1-theta[k,l])^(n[l]-x[i,l])
          }
          pz1[k] <- pi0[k]*prod(temp)
        }
        like[i] <- log(sum(pz1))
        like[like==-Inf] <- min(like[like!=-Inf])
      }
      likeli <- sum(like)
      likeli
    }

    ep <- 1
    num=0;
    logl <- NULL
    loglike0 <- 1
    while( ep > 1e-6){

      zik <- zik.new(x,theta,pi0,n,pathway)
      pi0 <- pi.new(zik)
      theta <- theta.new(x,zik,n)

      if (any(pi0<0.01)) {
        theta <- theta[pi0>=0.01,]
        pi0 <- pi0[pi0>=0.01]; pi0 <- pi0/sum(pi0)
      }

      (likeli <- llike(x,theta,pi0,n))

      (ep <- max(abs(likeli-loglike0)))

      loglike0 <- likeli
      num=num+1
      logl[num] <- likeli

    }

    C.hat <- apply(zik,1,which.max)
    prob <- apply(zik,1,max)
    BIC=-2*likeli+(length(pi0)+prod(dim(theta)))*log(dim(x)[1])
    list("z"=zik, "Pi"=pi0, "Theta"=theta,"Cluster"=C.hat, "Probability"=prob,"loglikelihood"=logl,BIC=BIC)
    }


  conv <- function(X,K,L,KCi,LCi,pathway=NULL){
    Y <- X
    N=dim(Y)[1];  p=dim(Y)[2];
    K=K; L=L;
    LC <- LCi
    KC <- KCi
	
    theta <- matrix(,K,L)
    for (k in 1:K){
      for (l in 1:L){
        theta[k,l] <-mean(Y[KC==k,LC==l])
      }
    }

    pi0 <- as.numeric(table(KC)/length(KC))

    x <- matrix(,dim(Y)[1],max(LC)); n <- NULL
    for (k in 1:max(LC)) {
      x[,k] <- apply(as.matrix(Y[,LC==k]),1,sum)
      n[k] <- sum(LC==k)
    }


    result <- list()
    num=0; ep=F; K1 <- KC; K2 <- LC

    while (ep==F) {
      num=num+1

      result[[1]] <- em(x,theta,pi0,n,pathway)

      Y <- t(Y)
      KC <- LC
      K <- max(KC)
      LC <- as.numeric(as.factor(result[[1]]$Cluster))
      L <- max(LC)

      x <- matrix(,dim(Y)[1],max(LC)); n <- NULL
      for (k in 1:max(LC)) {
        x[,k] <- apply(as.matrix(Y[,LC==k]),1,sum)
        n[k] <- sum(LC==k)
      }

      theta <- t(result[[1]]$Theta[sort(unique(result[[1]]$Cluster)),])
      pi0 <- as.numeric(table(KC)/length(KC))

      result[[2]] <- em(x,theta,pi0,n,NULL)

      Y <- t(Y)
      KC <- LC
      K <- max(KC)
      LC <- as.numeric(as.factor(result[[2]]$Cluster))
      L <- max(LC)

      x <- matrix(,dim(Y)[1],max(LC)); n <- NULL
      for (k in 1:max(LC)) {
        x[,k] <- apply(as.matrix(Y[,LC==k]),1,sum)
        n[k] <- sum(LC==k)
      }

      theta <- t(result[[2]]$Theta[sort(unique(result[[2]]$Cluster)),])
      pi0 <- as.numeric(table(KC)/length(KC))


      th1 <- result[[1]]$Theta
      th2 <- t(result[[2]]$Theta)

      ep1 <- all(K1==result[[1]]$Cluster)
      ep2 <- all(K2==result[[2]]$Cluster)
      K1 <- result[[1]]$Cluster
      K2 <- result[[2]]$Cluster
      ep <- all(ep1,ep2)

      if (num==20) break
    }
    result
  }

  Y <- X; N=dim(Y)[1]; p=dim(Y)[2]; B <- B; K=K; L=L

  if (is.null(path)) {
    KCi <- cutree(hclust(dist(Y),method="ward.D"),K)
    path.r <- path
  } else {

    path.kindex <- geneset <- NULL
    for (i in 1:length(path)) {
      path.k <- path[[i]]
      path.kindex <- c(path.kindex,match(path.k,rownames(Y)))
      geneset <- c(geneset,rep(i,length(path.k)))
    }

    Y <- rbind(Y[path.kindex,],Y[-path.kindex,])

    if (K==length(path)){
      path.r <- rep(K+1,N)
      path.r[1:length(path.kindex)] <- geneset

      KCi <- rep(K,N)
      KCi[1:length(path.kindex)] <- geneset
      geneset <- path.r

    } else {
      path.r <- rep(K,N)
      path.r[1:length(path.kindex)] <- geneset
      KCi <- path.r
      geneset <- path.r
    }

  }
  names(geneset) <- rownames(Y)

  LCi <- cutree(hclust(dist(t(Y)),method="ward.D"),L)

  result <- conv(Y,K,L,KCi,LCi,path.r)

  KC_candi <- result[[1]]$Cluster

  LC_candi <- result[[2]]$Cluster

  go_probability <- result[[2]]$Probability

  K1 <- length(unique(KC_candi))
  L1 <- length(unique(LC_candi))

  ### GO index for gene bootstrap
  geneindex <- goindex <- list()
  for (i in 1:B) {
    temp <- NULL
    for (j in 1:L1) {
      temp <- c(temp,sample(which(LC_candi==j),sum(LC_candi==j),replace=T))
    }
    goindex[[i]] <-  sort(temp)
  }

  message("---------------------\n")
  message("Starting bootstraps... \n")
  message("---------------------\n")

  genec1 <- matrix(,B,N)
  bKC <- KC_candi
  for (i in 1:B) {
    bdata <- Y[,goindex[[i]]]
    bLC <- LC_candi[goindex[[i]]]

    theta <- matrix(,K1,L1)
    for (k in 1:K1){
      for (l in 1:L1){
        theta[k,l] <-mean(bdata[bKC==k,bLC==l])
      }
    }
    pi0 <- as.numeric(table(bKC)/length(bKC))

    x1 <- matrix(,dim(bdata)[1],max(bLC)); n <- NULL
    for (k in 1:max(bLC)) {
      x1[,k] <- apply(as.matrix(bdata[,bLC==k]),1,sum)
      n[k] <- sum(bLC==k)
    }
    result <- em(x1,theta,pi0,n,path.r)
    genec1[i,] <- result$Cluster

    message(paste0("progress: ",i/B*100,"%"))
    flush.console()
  }

  temp1 <- apply(genec1,2,function(x) table(factor(x,levels=1:K)))/B

  KC <- apply(temp1,2,function(x) which.max(x))
  gene_pr <- apply(temp1,2,function(x) max(x))

  LC <- LC_candi
  go_pr <- go_probability

  names(KC) <- names(gene_pr) <- rownames(Y)
  names(LC) <- names(go_pr) <- colnames(Y)

  KC <- KC[match(rownames(X),names(KC))]
  gene_pr <- gene_pr[match(rownames(X),names(KC))]

  LC <- LC[match(colnames(X),names(LC))]
  go_pr <- go_pr[match(colnames(X),names(LC))]

  geneset <- geneset[match(rownames(X),names(geneset))]

  methods::new("palmer",
               data = X,
               init = list(nsample=dim(X)[1], nvariable=dim(X)[2], K = K, L = L, B = B),
               result =   list(Genecluster=KC, Geneprob=gene_pr, GOcluster=LC, GOprob=go_pr, Geneset=geneset)
  )
}

