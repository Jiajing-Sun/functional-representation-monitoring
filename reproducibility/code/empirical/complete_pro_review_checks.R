# Run from the manuscript directory. Reuses saved scores; no raw data changed.
source("reproducibility/code/R/dual_channel_monitor.R")
source("reproducibility/code/R/lrv_diagnostics.R")
dat <- function(f) read.csv(file.path("reproducibility/data",f),stringsAsFactors=FALSE)
a <- readRDS("reproducibility/code/empirical/outputs/spx_score_only_range_vs_shao_analysis.rds")
out <- "reproducibility/code/empirical/outputs/pro_review_checks"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
C <- a$fit$adjusted_range$prewhitener$C
ref <- a$fit$adjusted_range$calibration_sn_scores %*% t(C)
mon <- a$monitor$adjusted_range$sn_scores %*% t(C)
refdates <- as.Date(dat("spx_robustness_checks/spx_calibration_scores_2016_2019.csv")$date)
mondates <- as.Date(a$split_dates$monitor_dates)
cv <- dat("critical_value_grid_and_monte_carlo_checks.csv")
primary <- subset(cv,configuration=="exact997_251_empirical" & K==4)
one_alarm <- function(path,cut) { k<-which(path>cut); if(length(k)) as.character(mondates[k[1]]) else "None" }
uncertainty <- do.call(rbind,lapply(1:nrow(primary),function(i){
  p<-if(primary$method_id[i]=="adjusted_range") a$monitor$path$adjusted_range else a$monitor$path$shao
  data.frame(method=primary$method_id[i],critical=primary$critical_value[i],
    low=primary$quantile_ci95_low[i],high=primary$quantile_ci95_high[i],
    alarm_low=one_alarm(p,primary$quantile_ci95_low[i]),
    alarm_point=one_alarm(p,primary$critical_value[i]),
    alarm_high=one_alarm(p,primary$quantile_ci95_high[i]))
}))
write.csv(uncertainty,"reproducibility/data/spx_critical_value_alarm_sensitivity.csv",row.names=FALSE)
phi<-a$fit$adjusted_range$fpca$eigenfunctions[,1:4]
A<-solve(C); g<-phi%*%t(A)
write.csv(data.frame(u=(1:78-.5)/78,phi,g),"reproducibility/data/spx_original_and_demixed_loadings.csv",row.names=FALSE)
write.csv(A,"reproducibility/data/spx_inverse_ldl_matrix.csv",row.names=FALSE)
source("reproducibility/code/empirical/plot_spx_audit_figures.R")
plot_spx_demixed_loadings((1:78-.5)/78, phi, g,
  "figures/revised/spx_original_and_demixed_loadings.png")
# Repeat precisely the original preprocessing critical-value seeds and cache order.
pp<-dat("spx_robustness_checks/spx_preprocessing_sensitivity.csv")
keys<-paste(pp$K,pp$n_calibration,pp$n_monitor,sep="_")
first<-which(!duplicated(keys))
pcv<-parallel::mclapply(first,function(i){
  K<-pp$K[i];m<-pp$n_calibration[i];h<-pp$n_monitor[i];seed<-20260831L+i*10000L
  if(K==4 && m==997 && h==251) cuts<-primary$critical_value[match(c("adjusted_range","full_matrix"),primary$method_id)]
  else cuts<-c(dc_simulate_adjusted_range_critical_value(K,T=h/m,alpha=.05,nrep=20000,train_grid=m,seed=seed)$critical_value,
               dc_simulate_sn_critical_value(K,T=h/m,alpha=.05,nrep=20000,train_grid=m,seed=seed+500000L,ridge=0)$critical_value)
  data.frame(key=keys[i],range_cv=cuts[1],quadratic_cv=cuts[2],draws=if(K==4&&m==997)250000 else 20000)
},mc.cores=6)
pcv<-do.call(rbind,pcv);pp<-cbind(pp,pcv[match(keys,pcv$key),-1])
write.csv(pp,"reproducibility/data/spx_preprocessing_complete.csv",row.names=FALSE)
message("Recovered all preprocessing critical values.")
# Three historical windows specified before computing their paths. No window selection.
windows<-data.frame(name=c("2016--2018","2017--2019","2018--2019"),
                    from=c("2016-01-01","2017-01-01","2018-01-01"),
                    to=c("2018-12-31","2019-12-31","2019-12-31"))
draw_maxima <- function(job,m,h,K=4L,N=25000L) {
  set.seed(40400000L+job); s<-(1:h)/m; t<-(1:m)/m; ans<-matrix(0,N,2)
  for(i in 1:N){
    e<-matrix(rnorm((m+h)*K),m+h,K)/sqrt(m)
    B<-apply(e,2,cumsum); bridge<-B[1:m,,drop=FALSE]-outer(t,B[m,])
    ranges<-apply(rbind(0,bridge),2,function(x)diff(range(x)))
    U<-B[m+1:h,,drop=FALSE]-outer(1+s,B[m,])
    ans[i,1]<-max(rowSums(sweep(U,2,ranges,"/")^2)/(1+s)^2)
    vi<-chol2inv(chol(crossprod(bridge)/m))
    ans[i,2]<-max(rowSums((U%*%vi)*U)/(1+s)^2)
  }; ans
}
results<-list(); savedpaths<-list()
for(i in 1:nrow(windows)){
  sel<-refdates>=as.Date(windows$from[i]) & refdates<=as.Date(windows$to[i])
  rr<-ref[sel,,drop=FALSE];m<-nrow(rr);h<-nrow(mon)
  ci<-dc_ldl_factor(cov(rr)*(m-1)/m,ridge=0)$C
  rr<-rr%*%solve(t(ci));zz<-mon%*%solve(t(ci))
  paths<-list(adjusted_range=dc_sn_monitor_path(dc_fit_adjusted_range_sn(rr),zz),
              quadratic=dc_sn_monitor_path(dc_fit_full_matrix_sn(rr,ridge=0),zz))
  cache<-file.path(out,paste0("brownian_",m,"_",h,".rds"))
  if(file.exists(cache)) maxima<-readRDS(cache) else {
    maxima<-do.call(rbind,parallel::mclapply(i*100+1:10,draw_maxima,m=m,h=h,mc.cores=6))
    saveRDS(maxima,cache)
  }
  stopifnot(nrow(maxima)==250000,all(is.finite(maxima)))
  for(j in 1:2){
    cut<-quantile(maxima[,j],.95,type=8,names=FALSE)
    qq<-sort(maxima[,j]);ciidx<-pmin(length(qq),pmax(1,qbinom(c(.025,.975),length(qq),.95)+c(0,1)))
    H<-det_bartlett_lrv(rr,10)$lrv;corH<-cov2cor(H);diag(corH)<-0
    results[[length(results)+1]]<-data.frame(reference=windows$name[i],m=m,h=h,K=4,
      method=names(paths)[j],critical=cut,low=qq[ciidx[1]],high=qq[ciidx[2]],
      first_alarm=one_alarm(paths[[j]],cut),maximum_ratio=max(paths[[j]]/cut),
      offdiag_diagnostic=sqrt(sum(corH^2)/4),draws=nrow(maxima))
    savedpaths[[length(savedpaths)+1]]<-data.frame(reference=windows$name[i],method=names(paths)[j],
                                                date=mondates,statistic=paths[[j]],critical=cut)
  }
  message("Completed historical reference window ",windows$name[i])
}
write.csv(do.call(rbind,results),"reproducibility/data/spx_reference_window_sensitivity.csv",row.names=FALSE)
write.csv(do.call(rbind,savedpaths),"reproducibility/data/spx_reference_window_paths.csv",row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(uncertainty)
print(do.call(rbind,results))
