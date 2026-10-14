# Run from the manuscript directory. Prespecified supplementary null experiment.
# The fit regime differs from the stationary post-fit regime; no DGP search.
source("reproducibility/code/R/dual_channel_monitor.R")
source("reproducibility/code/simulations/formal/R/null_dgp.R")
args <- commandArgs(TRUE)
nrep <- if (length(args)) as.integer(args[1]) else 5000L
workers <- if (length(args)>1) as.integer(args[2]) else 6L
outdir <- "reproducibility/code/simulations/componentwise_range/outputs/nonidentity_ldl"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
Q <- dc_functional_basis_rotation(8L)
write.csv(Q, "reproducibility/data/fixed_functional_rotation_8.csv", row.names=FALSE)
cv <- read.csv("reproducibility/data/critical_value_grid_and_monte_carlo_checks.csv")
cv <- subset(cv, configuration=="exact1000_500_T0.5")
design <- data.frame(id=1:6, n=c(500,500,500,250,1000,500),
                     K=c(2,2,2,2,2,4), mixing=c(0,.5,1,1,1,1))
methods <- c("no_LDL", "fitted_LDL", "oracle_LDL", "quadratic")
one <- function(rep, cfg) {
  set.seed(20260905L + 100000L*cfg$id + rep)
  n<-cfg$n; K<-cfg$K; m<-1000L; h<-500L; N<-m+h
  C<-diag(8); C[cbind(2:8,1:7)]<-cfg$mixing
  lambda<-dc_score_eigenvalues(8)
  xf<-sweep(matrix(rnorm(n*8),n,8),2,sqrt(lambda),"*") %*% t(Q)
  ef<-eigen(cov(xf),symmetric=TRUE)$vectors[,1:K,drop=FALSE]
  # Sign alignment preserves the prespecified population FPCA order.
  sg<-sign(diag(crossprod(Q[,1:K,drop=FALSE],ef))); sg[sg==0]<-1
  ef<-sweep(ef,2,sg,"*")
  ph<-seq(.15,.55,length.out=8)
  u<-matrix(0,N,8)
  for(j in 1:8) u[,j]<-as.numeric(filter(rnorm(N,sd=sqrt(1-ph[j]^2)),
                                         ph[j],method="recursive",init=rnorm(1)))
  scores<-u %*% t(C) %*% t(Q) %*% ef
  ref<-scores[1:m,,drop=FALSE]; mon<-scores[m+1:h,,drop=FALSE]
  A<-dc_ldl_factor(cov(ref)*(m-1)/m,ridge=0)$C
  fitted_ref<-ref %*% solve(t(A)); fitted_mon<-mon %*% solve(t(A))
  # The oracle knows both the limiting fit eigenfunctions and population LDL.
  oracle_ref<-u[1:m,1:K,drop=FALSE]; oracle_mon<-u[m+1:h,1:K,drop=FALSE]
  paths<-list(
    no_LDL=dc_sn_monitor_path(dc_fit_adjusted_range_sn(ref),mon),
    fitted_LDL=dc_sn_monitor_path(dc_fit_adjusted_range_sn(fitted_ref),fitted_mon),
    oracle_LDL=dc_sn_monitor_path(dc_fit_adjusted_range_sn(oracle_ref),oracle_mon),
    quadratic=dc_sn_monitor_path(dc_fit_full_matrix_sn(fitted_ref,ridge=0),fitted_mon))
  c(rep=rep, vapply(paths,max,numeric(1)))
}
summaries<-list(); population<-list()
for(ii in 1:nrow(design)) {
  cfg<-design[ii,]
  rows<-parallel::mclapply(1:nrep,one,cfg=cfg,mc.cores=workers)
  stopifnot(all(vapply(rows,is.numeric,logical(1))))
  raw<-as.data.frame(do.call(rbind,rows)); stopifnot(all(is.finite(as.matrix(raw))))
  saveRDS(list(design=cfg,nrep=nrep,raw=raw),file.path(outdir,paste0("design_",ii,".rds")))
  cuts<-setNames(c(rep(cv$critical_value[cv$K==cfg$K & cv$method_id=="adjusted_range"],3),
                  cv$critical_value[cv$K==cfg$K & cv$method_id=="full_matrix"]),methods)
  for(nm in methods) {
    y<-raw[[nm]]>cuts[nm]; z<-raw$fitted_LDL>cuts["fitted_LDL"]
    summaries[[length(summaries)+1]]<-data.frame(cfg,m=1000,h=500,nrep=nrep,method=nm,
      critical=unname(cuts[nm]),rejection=mean(y),mcse=sqrt(mean(y)*(1-mean(y))/nrep),
      difference_from_fitted=mean(y-z),paired_se=sd(as.integer(y)-as.integer(z))/sqrt(nrep))
  }
  K<-cfg$K; C<-diag(K); if(K>1) C[cbind(2:K,1:(K-1))]<-cfg$mixing
  ph<-seq(.15,.55,length.out=8)[1:K]; Om<-diag((1+ph)/(1-ph))
  Gamma<-C%*%t(C); Sigma<-C%*%Om%*%t(C)
  population[[ii]]<-list(design=cfg,C=C,Gamma=Gamma,Sigma=Sigma,Omega=Om,
                         raw_LRV_correlation=cov2cor(Sigma),post_LDL_correlation=diag(K))
  message("Completed design ",ii,"; n=",cfg$n,", K=",K,", mixing=",cfg$mixing)
}
write.csv(do.call(rbind,summaries),"reproducibility/data/nonidentity_ldl_null_5000rep.csv",row.names=FALSE)
saveRDS(population,file.path(outdir,"population_matrices.rds"))
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
