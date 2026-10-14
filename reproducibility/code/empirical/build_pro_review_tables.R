# Run from the manuscript directory after complete_pro_review_checks.R.
dat<-function(f)read.csv(file.path("reproducibility/data",f),stringsAsFactors=FALSE)
fmt<-function(x,n=4)sprintf(paste0("%.",n,"f"),x)
row<-function(...)paste0(paste(...,collapse=" & ")," \\\\")
tabletex<-function(file,label,caption,cols,header,body,note){
  writeLines(c("\\begin{table}[!htbp]","\\centering\\small\\color{revisiondeepred}",
    "\\setlength{\\tabcolsep}{4pt}",paste0("\\caption{",caption,"}\\label{",label,"}"),
    paste0("\\begin{tabular}{",cols,"}"),"\\toprule",header,"\\midrule",body,
    "\\bottomrule\\end{tabular}",paste0("\\begin{flushleft}\\small Note: ",note,"\\end{flushleft}"),
    "\\end{table}"),file.path("tables",file))
}
cv<-dat("critical_value_grid_and_monte_carlo_checks.csv")
d<-subset(cv,configuration=="exact1000_500_T0.5")
rows<-sapply(1:5,function(k){r<-subset(d,K==k&method_id=="adjusted_range");q<-subset(d,K==k&method_id=="full_matrix")
 row(c(k,paste0("[",fmt(r$quantile_ci95_low),", ",fmt(r$quantile_ci95_high),"]"),
          paste0("[",fmt(q$quantile_ci95_low),", ",fmt(q$quantile_ci95_high),"]")))})
tabletex("pro_critical_intervals.tex","tab:critical-intervals","Monte Carlo uncertainty of matched-grid critical values.","ccc",
 row(c("\\(K\\)","Adjusted-range SN, 95\\% interval","Quadratic SN, 95\\% interval")),rows,
 "The reference and monitoring grids have 1,000 and 500 points. Intervals are binomial order-statistic intervals from 250,000 draws, not confidence intervals for the data-generating process.")
ldl<-dat("nonidentity_ldl_null_5000rep.csv")
rows<-sapply(unique(ldl$id),function(id){
 d<-ldl[ldl$id==id,];row(c(d$n[1],d$K[1],fmt(d$mixing[1],1),fmt(d$rejection[match(c("no_LDL","fitted_LDL","oracle_LDL","quadratic"),d$method)])))})
tabletex("pro_nonidentity_ldl.tex","tab:nonidentity-ldl","Null rejection when the fit and post-fit covariances differ.","rrrcccc",
 row(c("\\(n\\)","\\(K\\)","\\(c\\)","No LDL","Fitted LDL","Oracle LDL","Quadratic SN")),rows,
 "Each row uses 5,000 paired replications, with reference length 1,000, horizon 500 and nominal level 5\\%. The first three columns of results use adjusted-range SN. The oracle knows both the limiting fit basis and the population LDL map.")
ann<-dat("spx_robustness_checks/spx_calibration_annual_summary.csv")
cal<-dat("spx_robustness_checks/spx_calibration_scores_2016_2019.csv")
sc<-apply(cal[,-1],2,sd)
rows<-sapply(2016:2019,function(y){d<-ann[ann$year==y,];d<-d[order(d$score),]
 row(c(y,d$n[1],paste0(fmt(d$mean/sc,3)," (",fmt(d$sd/sc,3),")")))})
tabletex("pro_reference_annual.tex","tab:reference-annual","Annual summaries of all four reference scores.","crcccc",
 row(c("Year","\\(n\\)","Score 1","Score 2","Score 3","Score 4")),rows,
 "Each entry is the annual mean (standard deviation), divided by that coordinate's full-reference standard deviation. Scores are in the original FPCA coordinates and centered over 2016--2019.")
pp<-dat("spx_preprocessing_complete.csv")
datefmt<-function(x)ifelse(is.na(x)|!nzchar(x),"None",substr(x,6,10))
for(part in 1:2){
 d<-pp[((part-1)*10+1):(part*10),]
 rows<-sapply(1:nrow(d),function(i){x<-d[i,];row(c(fmt(100*x$fve,0),x$spline_df,x$bin_minutes,
  if(x$min_day_coverage==1)"100" else "95",x$K,x$n_calibration,
  datefmt(x$adjusted_range_alarm),datefmt(x$full_matrix_alarm),
  fmt(x$adjusted_range_max_ratio,3),fmt(x$full_matrix_max_ratio,3)))})
 tabletex(paste0("pro_preprocessing_",part,".tex"),paste0("tab:preprocessing-",part),
   paste0("Complete preprocessing comparison (",part," of 2)."),"rrrrrrllrr",
   row(c("FVE","Spline","Bin","Coverage","\\(K\\)","\\(m\\)","Range date","Quad. date","Range max","Quad. max")),rows,
   "FVE and coverage are percentages; bin width is in minutes. All rows have fit length 497 and monitoring horizon 251. Dates are month--day in 2020; None means no crossing by the horizon. Maximum entries divide each statistic by its own critical value. Critical values are in Table~\\ref{tab:preprocessing-cv}.")
}
d<-unique(pp[,c("K","n_calibration","n_monitor","range_cv","quadratic_cv","draws")]);d<-d[order(d$n_calibration,d$K),]
rows<-sapply(1:nrow(d),function(i){x<-d[i,];row(c(x$K,x$n_calibration,x$n_monitor,fmt(x$range_cv),fmt(x$quadratic_cv),format(x$draws,big.mark=",")))})
tabletex("pro_preprocessing_cv.tex","tab:preprocessing-cv","Critical values used in the preprocessing comparison.","rrrccr",
 row(c("\\(K\\)","\\(m\\)","\\(h\\)","Adjusted-range SN","Quadratic SN","Draws")),rows,
 "The original simulation seeds and cache order have been retained. The earlier 2020-only configuration uses 250,000 draws; auxiliary configurations use 20,000. Each row uses its own matched finite grid.")
bs<-dat("spx_fit_block_bootstrap_stability.csv");bs<-subset(bs,success)
rows<-sapply(sort(unique(bs$K_selected)),function(k){d<-bs[bs$K_selected==k,]
 r<-!is.na(d$adjusted_range_alarm_index);q<-!is.na(d$quadratic_alarm_index);both<-r&q
 row(c(k,nrow(d),sum(r),sum(q),sum(both),sum(d$adjusted_range_alarm_index[both]<d$quadratic_alarm_index[both])))})
tabletex("pro_bootstrap.tex","tab:bootstrap-complete","Fit-block bootstrap counts by selected dimension.","rrrrrr",
 row(c("\\(\\widehat K\\)","Fits","Range alarms","Quad. alarms","Both","Range earlier")),rows,
 "There are 999 successful fits and one failure to reach 90\\% FVE at the cap of five. The last column counts strict earlier crossings among paths on which both procedures signal. These are conditional sensitivity counts, not probabilities of a true structural change.")
d<-dat("spx_reference_window_sensitivity.csv")
rows<-sapply(unique(d$reference),function(w){x<-d[d$reference==w,];x<-x[match(c("adjusted_range","quadratic"),x$method),]
 row(c(w,x$m[1],fmt(x$critical),datefmt(x$first_alarm),fmt(x$maximum_ratio,3)))})
tabletex("pro_reference_windows.tex","tab:reference-windows","2020 monitoring with alternative historical reference windows.","lrcc ll cc",
 row(c("Reference","\\(m\\)","Range CV","Quad. CV","Range date","Quad. date","Range max","Quad. max")),rows,
 "All rows keep the 2014--2015 basis, four coordinates and 251-day 2020 horizon fixed. The center, LDL map and both normalizers are re-estimated within each reference window. Every critical value uses 250,000 matched-grid draws. These earlier comparisons monitor 2020 alone; they are separate from the continuous 2019--2020 analysis.")
d<-dat("spx_critical_value_alarm_sensitivity.csv")
rows<-sapply(1:nrow(d),function(i){x<-d[i,];row(c(if(i==1)"Adjusted-range SN" else "Quadratic SN",
  paste0("[",fmt(x$low),", ",fmt(x$high),"]"),datefmt(x$alarm_low),datefmt(x$alarm_point),datefmt(x$alarm_high)))})
tabletex("pro_alarm_uncertainty.tex","tab:alarm-uncertainty","Alarm dates at Monte Carlo interval endpoints.","lclll",
 row(c("Method","Critical-value 95\\% interval","Lower","Point","Upper")),rows,
 "Dates correspond to substituting the lower endpoint, point estimate and upper endpoint into the same earlier 2020-only path. This measures only critical-value simulation error; it is not an alarm-date confidence interval.")

p<-subset(dat("paired_power_summary_5000rep.csv"),calibration=="size_adjusted")
diffs<-subset(dat("paired_power_differences_5000rep.csv"),calibration=="size_adjusted")
dgps<-c(gaussian_iid="DGP1",hong_var_diag="DGP2",component_garch="DGP4",component_sv="DGP5",component_t5="DGP6")
rows<-sapply(1:nrow(diffs),function(i){x<-diffs[i,];s<-p[p$scenario_id==x$scenario_id,]
 rr<-s$prechange_false_alarm[match(c("adjusted_range","shao"),s$method_id)]
 # The archived benchmark identifier is full_matrix in some versions.
 if(anyNA(rr)) rr<-s$prechange_false_alarm[match(c("adjusted_range","full_matrix"),s$method_id)]
 stopifnot(length(rr)==2,all(is.finite(rr)))
 row(c(dgps[x$design],if(x$change_shape=="abrupt")"Abrupt" else "Ramp",
       if(x$direction=="sparse_first")"Sparse" else "Dense",fmt(x$delta_standardized,2),
       fmt(rr),x$common_risk_n))})
header<-row(c("Design","Change","Direction","Signal","Range pre-alarm","Quad. pre-alarm","Common risk"))
writeLines(c("{\\small\\color{revisiondeepred}\\setstretch{1.05}",
 "\\setlength{\\tabcolsep}{4pt}\\begin{longtable}{lllrrrr}",
 "\\caption{Pre-change alarms and common risk-set counts for the size-adjusted power comparison.}\\label{tab:power-risk}\\\\",
 "\\toprule",header,"\\midrule\\endfirsthead","\\toprule",header,"\\midrule\\endhead",rows,
 "\\bottomrule\\end{longtable}",
 "\\noindent Each cell contains 5,000 replications. Pre-alarm is the probability of crossing by monitoring observation 100. Common risk counts replications with no such crossing under either method; capped mean delays use this shared set.}",
 ""),"tables/pro_power_risk_set.tex")
