#!/usr/bin/env Rscript
# Run from the manuscript directory. Fixed design; no window selection.
# Reuse only the 2014--2015 FPCA fit. Re-estimate every reference quantity.
source("reproducibility/code/R/dual_channel_monitor.R")
source("reproducibility/code/R/lrv_diagnostics.R")
out <- "reproducibility/code/empirical/outputs/continuous_2019_2020"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
save_csv <- function(x, suffix) write.csv(x, paste0("reproducibility/data/spx_continuous_",suffix,".csv"),row.names=FALSE,na="")
a <- readRDS("reproducibility/code/empirical/outputs/spx_score_only_range_vs_shao_analysis.rds")
old <- a$fit$adjusted_range
stopifnot(old$fpca$K==4L, old$feature_prewhitener$order==0L)
# Undo the old common affine transformation exactly. Its old reference center
# cancels on recentering in the new reference block (verified below).
recover <- function(z) sweep(z %*% t(old$prewhitener$C),2,old$prewhitener$center,"+")
scores <- rbind(recover(old$calibration_sn_scores),recover(a$monitor$adjusted_range$sn_scores))
dates <- as.Date(c(a$split_dates$calibration_dates,a$split_dates$monitor_dates))
stopifnot(!anyDuplicated(dates),all(diff(dates)>0),nrow(scores)==1248L)
colnames(scores)<-paste0("PC",1:4)
save_csv(data.frame(date=dates,scores),"scores")
select <- function(from,to) which(dates>=as.Date(from)&dates<=as.Date(to))
fit_paths <- function(ref_id,mon_id) {
  rr<-scores[ref_id,,drop=FALSE];xx<-scores[mon_id,,drop=FALSE]
  ldl<-dc_fit_ldl_prewhitener(rr,ridge=0,mode="demix")
  z<-dc_apply_ldl_prewhitener(ldl,rr);w<-dc_apply_ldl_prewhitener(ldl,xx)
  sn<-list(adjusted_range=dc_fit_adjusted_range_sn(z,ridge=0),quadratic=dc_fit_full_matrix_sn(z,ridge=0))
  paths<-vapply(sn,dc_sn_monitor_path,numeric(nrow(w)),monitor_features=w)
  # Independent expression for the statistic and invariance to a common offset.
  U<-apply(sweep(w,2,colMeans(z),"-"),2,cumsum)/sqrt(nrow(z))
  boundary<-(1+seq_len(nrow(w))/nrow(z))^2
  direct<-cbind(rowSums(sweep(U,2,sn$adjusted_range$adjusted_ranges,"/")^2)/boundary,
                rowSums((U%*%solve(sn$quadratic$self_normalizer))*U)/boundary)
  stopifnot(max(abs(paths-direct))<1e-8,ldl$ridge_added==0,
            sn$quadratic$ridge_added==0,all(sn$adjusted_range$adjusted_ranges>1e-7))
  qraw<-dc_sn_monitor_path(dc_fit_full_matrix_sn(rr,ridge=0),xx)
  stopifnot(max(abs(qraw-paths[,2]))<1e-8)
  list(ref=rr,monitor=xx,z=z,w=w,ldl=ldl,sn=sn,paths=paths,U=U,
       dates=dates[mon_id],ref_dates=dates[ref_id],m=nrow(rr),h=nrow(xx))
}
designs<-list(
  continuous=fit_paths(select("2016-01-01","2018-12-31"),select("2019-01-01","2020-12-31")),
  reference_2018=fit_paths(select("2016-01-01","2017-12-31"),select("2018-01-01","2018-12-31")),
  reference_2017_2018=fit_paths(select("2016-01-01","2016-12-31"),select("2017-01-01","2018-12-31")))
stopifnot(designs$continuous$m==748,designs$continuous$h==500)
# The 2019 prefix must equal the previously saved one-year statistic before
# division by a critical value; only the longer-horizon threshold changes.
legacy<-read.csv("reproducibility/data/spx_robustness_checks/spx_2019_placebo_path.csv")
oldcv<-subset(read.csv("reproducibility/data/critical_value_grid_and_monte_carlo_checks.csv"),
               configuration=="exact748_249_placebo" & K==4)
legacycuts<-oldcv$critical_value[match(c("adjusted_range","full_matrix"),oldcv$method_id)]
stopifnot(nrow(legacy)==249,length(legacycuts)==2,
  max(abs(designs$continuous$paths[1:249,1]/legacycuts[1]-legacy$adjusted_range_ratio))<1e-8,
  max(abs(designs$continuous$paths[1:249,2]/legacycuts[2]-legacy$full_matrix_ratio))<1e-8)
saveRDS(designs,file.path(out,"fitted_designs.rds"))

# Exact Brownian values at integer/m grid nodes. Same draws for both methods;
# reference bridge includes r=0 for the ranges and right endpoints for V.
draw_chunk<-function(job,m,h,seed,N=5000L,prefix=0L) {
  set.seed(seed+job);K<-4L;s<-(1:h)/m;t<-(1:m)/m
  ans<-matrix(NA_real_,N,if(prefix>0)4L else 2L)
  for(i in seq_len(N)) {
    B<-apply(matrix(rnorm((m+h)*K),m+h,K)/sqrt(m),2,cumsum)
    bridge<-B[1:m,,drop=FALSE]-outer(t,B[m,])
    ranges<-apply(rbind(0,bridge),2,function(x)diff(range(x)))
    U<-B[m+seq_len(h),,drop=FALSE]-outer(1+s,B[m,])
    r<-rowSums(sweep(U,2,ranges,"/")^2)/(1+s)^2
    inv<-chol2inv(chol(crossprod(bridge)/m))
    q<-rowSums((U%*%inv)*U)/(1+s)^2
    ans[i,1:2]<-c(max(r),max(q))
    if(prefix>0) ans[i,3:4]<-c(max(r[1:prefix]),max(q[1:prefix]))
  }
  ans
}
cores<-6L
cv_rows<-list();alarm_rows<-list();all_paths<-list();contrib<-list();diag_rows<-list()
for(ii in seq_along(designs)) {
  name<-names(designs)[ii];d<-designs[[ii]];cache<-file.path(out,paste0(name,"_brownian_250000.rds"))
  seed<-202609070L+ii*1000L
  if(file.exists(cache)) maxima<-readRDS(cache) else {
    pieces<-list()
    for(batch in 0:4) {
      jobs<-batch*10+1:10
      pieces<-c(pieces,parallel::mclapply(jobs,draw_chunk,m=d$m,h=d$h,seed=seed,
                  prefix=if(name=="continuous")249L else 0L,mc.cores=cores,mc.set.seed=FALSE))
      message(name,": ",(batch+1)*50000," Brownian draws completed")
    }
    maxima<-do.call(rbind,pieces);saveRDS(maxima,cache)
  }
  stopifnot(nrow(maxima)==250000L,all(is.finite(maxima)))
  if(name=="continuous") stopifnot(all(maxima[,1:2]>=maxima[,3:4]))
  cuts<-apply(maxima[,1:2],2,quantile,probs=.95,type=8,names=FALSE)
  first<-function(path,cut) {k<-which(path>cut);if(length(k)) k[1] else NA_integer_}
  for(j in 1:2) {
    method<-colnames(d$paths)[j];q<-sort(maxima[,j]);ix<-pmin(length(q),pmax(1,qbinom(c(.025,.975),length(q),.95)+c(0,1)))
    lo<-q[ix[1]];hi<-q[ix[2]];path<-d$paths[,j];ratio<-path/cuts[j];k<-first(path,cuts[j])
    cv_rows[[length(cv_rows)+1]]<-data.frame(design=name,method=method,K=4,m=d$m,h=d$h,T=d$h/d$m,
        draws=nrow(maxima),alpha=.05,critical=cuts[j],ci_low=lo,ci_high=hi,seed_base=seed)
    alarm_rows[[length(alarm_rows)+1]]<-data.frame(design=name,method=method,alarm_index=k,
        alarm_date=as.character(d$dates[k]),ratio_at_alarm=ratio[k],maximum_ratio=max(ratio),
        maximum_date=as.character(d$dates[which.max(ratio)]),year_end_2019=ratio[match(as.Date("2019-12-31"),d$dates)],
        march16_2020=ratio[match(as.Date("2020-03-16"),d$dates)],final_ratio=tail(ratio,1),
        alarm_low_cv=as.character(d$dates[first(path,lo)]),alarm_high_cv=as.character(d$dates[first(path,hi)]))
    all_paths[[length(all_paths)+1]]<-data.frame(design=name,date=d$dates,k=seq_len(d$h),method=method,statistic=path,critical=cuts[j],ratio=ratio)
    if(name=="continuous" && !is.na(k)) {
      v<-d$U[k,];terms<-v*as.numeric(d$sn[[j]]$inverse%*%v)
      shift<-colMeans(d$monitor[seq_len(k),,drop=FALSE])-colMeans(d$ref)
      contrib[[j]]<-data.frame(method=method,date=d$dates[k],coordinate=1:4,original_score_mean_shift=shift,
          demixed_score_mean_shift=as.numeric(solve(d$ldl$C)%*%shift),fraction=terms/sum(terms),
          ratio_contribution=terms/(1+k/d$m)^2/cuts[j])
      stopifnot(abs(sum(terms)/(1+k/d$m)^2-path[k])<1e-8)
    }
  }
  for(bw in c(5,10,20,40)) {
    P<-cov2cor(det_bartlett_lrv(d$z,bw)$lrv);diag(P)<-0
    diag_rows[[length(diag_rows)+1]]<-data.frame(design=name,bandwidth=bw,diagnostic=sqrt(sum(P^2)/4),max_abs=max(abs(P)))
  }
  d$critical<-cuts;designs[[ii]]<-d
  message(name," complete; critical values ",paste(round(cuts,6),collapse=", "))
}
save_csv(do.call(rbind,cv_rows),"critical_values")
save_csv(do.call(rbind,alarm_rows),"alarms")
save_csv(do.call(rbind,all_paths),"paths")
save_csv(do.call(rbind,contrib),"contributions")
save_csv(do.call(rbind,diag_rows),"correlation")
d<-designs$continuous
annual<-do.call(rbind,lapply(c("2016","2017","2018"),function(year) {
  ii<-format(d$ref_dates,"%Y")==year
  data.frame(year=year,coordinate=1:4,n=sum(ii),mean=colMeans(d$ref[ii,,drop=FALSE])-colMeans(d$ref),
    sd=apply(d$ref[ii,,drop=FALSE],2,sd),full_reference_sd=apply(d$ref,2,sd))
}))
save_csv(annual,"reference_annual")
save_csv(data.frame(coordinate=1:4,acf1=vapply(1:4,function(j)acf(d$z[,j],lag.max=20,plot=FALSE)$acf[2],numeric(1)),
   acf20=vapply(1:4,function(j)acf(d$z[,j],lag.max=20,plot=FALSE)$acf[21],numeric(1))),"reference_acf")
save_csv(data.frame(method=names(d$sn),condition_number=vapply(d$sn,function(x)x$condition_number,numeric(1)),
  ridge_added=vapply(d$sn,function(x)x$ridge_added,numeric(1))),"normalizers")
save_csv(data.frame(coordinate=1:4,solve(d$ldl$C)),"inverse_ldl")
phi<-old$fpca$eigenfunctions[,1:4];g<-phi%*%t(solve(d$ldl$C))
save_csv(data.frame(u=(1:78-.5)/78,phi,g),"loadings")
saveRDS(designs,file.path(out,"complete_analysis.rds"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(do.call(rbind,cv_rows));print(do.call(rbind,alarm_rows));print(do.call(rbind,diag_rows))
