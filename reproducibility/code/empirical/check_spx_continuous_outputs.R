# Independent checks of the saved continuous-monitoring results.
source("reproducibility/code/R/dual_channel_monitor.R")
out<-"reproducibility/code/empirical/outputs/continuous_2019_2020"
designs<-readRDS(file.path(out,"complete_analysis.rds"))
cv<-read.csv("reproducibility/data/spx_continuous_critical_values.csv")
for(name in names(designs)) {
  d<-designs[[name]];c<-subset(cv,design==name)
  draws<-readRDS(file.path(out,paste0(name,"_brownian_250000.rds")))
  set.seed(c$seed_base[1]+1L)
  for(i in 1:10) {
    x<-matrix(rnorm((d$m+d$h)*4),d$m+d$h,4)
    rr<-x[1:d$m,,drop=FALSE];zz<-x[d$m+seq_len(d$h),,drop=FALSE]
    r<-dc_sn_monitor_path(dc_fit_adjusted_range_sn(rr),zz)
    q<-dc_sn_monitor_path(dc_fit_full_matrix_sn(rr),zz)
    stopifnot(max(abs(c(max(r),max(q))-draws[i,1:2]))<1e-9)
    if(name=="continuous")stopifnot(max(abs(c(max(r[1:249]),max(q[1:249]))-draws[i,3:4]))<1e-9)
  }
  # Reconstruct paths after adding a common arbitrary shift to all scores.
  offset<-c(10,-4,2,7)
  ref<-sweep(d$ref,2,offset,"+");mon<-sweep(d$monitor,2,offset,"+")
  ldl<-dc_fit_ldl_prewhitener(ref,ridge=0,mode="demix")
  rr<-dc_apply_ldl_prewhitener(ldl,ref);zz<-dc_apply_ldl_prewhitener(ldl,mon)
  p<-cbind(dc_sn_monitor_path(dc_fit_adjusted_range_sn(rr),zz),dc_sn_monitor_path(dc_fit_full_matrix_sn(rr),zz))
  stopifnot(max(abs(p-d$paths))<1e-8)
  message(name,": Brownian draws and centering invariance verified")
}
d<-designs$continuous
alarms<-subset(read.csv("reproducibility/data/spx_continuous_alarms.csv"),design=="continuous")
con<-read.csv("reproducibility/data/spx_continuous_contributions.csv")
phi<-readRDS("reproducibility/code/empirical/outputs/spx_score_only_range_vs_shao_analysis.rds")$fit$adjusted_range$fpca$eigenfunctions[,1:4]
for(j in 1:2){
  k<-alarms$alarm_index[j]
  stopifnot(all(d$paths[seq_len(k-1),j]<d$critical[j]),d$paths[k,j]>d$critical[j])
  shift<-colMeans(d$monitor[seq_len(k),,drop=FALSE])-colMeans(d$ref)
  message(names(d$sn)[j]," reconstructed log-RV contrast range: ",paste(range(phi%*%shift),collapse=", "))
  stopifnot(abs(sum(con$fraction[con$method==names(d$sn)[j]])-1)<1e-10)
}
message("Maximum 2019 ratios: ",paste(apply(d$paths[1:249,,drop=FALSE],2,max)/d$critical,collapse=", "))
message("All continuous-monitoring checks passed.")
