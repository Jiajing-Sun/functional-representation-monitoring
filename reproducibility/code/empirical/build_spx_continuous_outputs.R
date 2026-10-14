#!/usr/bin/env Rscript
# Presentation only. Run after run_spx_continuous_2019_2020.R.
read_result<-function(s) read.csv(paste0("reproducibility/data/spx_continuous_",s,".csv"),stringsAsFactors=FALSE)
cv<-read_result("critical_values");alarms<-read_result("alarms");paths<-read_result("paths")
names_method<-c(adjusted_range="Adjusted-range SN",quadratic="Quadratic SN")
date_text<-function(x) {
  d<-as.Date(x);month<-format(d,"%b");month<-ifelse(month=="May",month,paste0(month,"."))
  ifelse(is.na(d),"None",paste0(month,"~",as.integer(format(d,"%d")),", ",format(d,"%Y")))
}
tex_table<-function(file,label,caption,align,header,rows,note) {
  scale_table<-!grepl("checkpoints|annual",file)
  writeLines(c("\\begin{table}[!htbp]","\\centering\\small", "\\color{revisioncontinuousblue}",
    paste0("\\caption{\\contrev{",caption,"}}\\label{",label,"}"),"\\setlength{\\tabcolsep}{4pt}",
    if(scale_table)"\\resizebox{\\linewidth}{!}{%",paste0("\\begin{tabular}{",align,"}"),"\\toprule",paste0(header," \\\\"),
    "\\midrule",paste0(rows," \\\\"),paste0("\\bottomrule\\end{tabular}",if(scale_table)"}" else ""),
    paste0("\\begin{flushleft}\\small Note: ",note,"\\end{flushleft}"),"\\end{table}"),file)
}
cmain<-subset(cv,design=="continuous");amain<-subset(alarms,design=="continuous")
tex_table("tables/spx_continuous_main.tex","tab:spx","Continuous S\\&P~500 monitoring over 2019--2020.","lrrlrr",
  "Method & Critical value & Alarm index & First alarm & End-2019 ratio & End-2020 ratio",
  vapply(1:2,function(i)sprintf("%s & %.4f & %d & %s & %.4f & %.4f",names_method[cmain$method[i]],cmain$critical[i],
    amain$alarm_index[i],date_text(amain$alarm_date[i]),amain$year_end_2019[i],amain$final_ratio[i]),character(1)),
  "The reference sample is 2016--2018 ($m=748$); 500 observations are monitored without resetting at the start of 2020. Each ratio uses its method-specific 5\\% critical value for the entire two-year horizon. Values after a first crossing describe the unchanged statistic, not a new alarm test.")
design_names<-c(continuous="2019--2020",reference_2018="2018",reference_2017_2018="2017--2018")
tex_table("tables/spx_continuous_cv.tex","tab:continuous-cv","Matched-grid critical values for continuous monitoring and reference-period checks.","llrrrc",
  "Monitoring period & Method & $m$ & $h$ & Critical value & MC 95\\% interval",
  vapply(1:nrow(cv),function(i)sprintf("%s & %s & %d & %d & %.4f & [%.4f, %.4f]",design_names[cv$design[i]],
   names_method[cv$method[i]],cv$m[i],cv$h[i],cv$critical[i],cv$ci_low[i],cv$ci_high[i]),character(1)),
  "Every row uses $K=4$, a 5\\% nominal level and 250,000 matched-grid Brownian draws. Within a design, the methods use the same draws. The intervals quantify Monte Carlo error in the quantile, not uncertainty about reference stability.")
aref<-subset(alarms,design!="continuous")
tex_table("tables/spx_continuous_reference.tex","tab:continuous-reference","Sequential checks within the proposed 2016--2018 reference period.","llllr",
  "Reference & Monitored period & Method & First crossing & Maximum ratio",
  vapply(1:nrow(aref),function(i)sprintf("%s & %s & %s & %s & %.4f",
    if(aref$design[i]=="reference_2018")"2016--2017" else "2016",design_names[aref$design[i]],names_method[aref$method[i]],
    date_text(aref$alarm_date[i]),aref$maximum_ratio[i]),character(1)),
  "Both checks keep the 2014--2015 FPCA fit and $K=4$, re-estimating reference quantities before each monitored period. Each method uses its own 5\\% critical value for that complete horizon. These overlapping checks are reported separately, not as independent tests or a joint 5\\% decision.")
con<-read_result("contributions")
tex_table("tables/spx_continuous_contributions.tex","tab:continuous-contributions","Score-mean changes and statistic contributions at the continuous-monitoring alarms.","llrrr",
  "Method & Coordinate & Original-score mean change & Demixed-score mean change & Contribution (\\%)",
  vapply(1:nrow(con),function(i)sprintf("%s & %d & %.4f & %.4f & %.2f",names_method[con$method[i]],con$coordinate[i],
    con$original_score_mean_shift[i],con$demixed_score_mean_shift[i],100*con$fraction[i]),character(1)),
  "Mean changes compare cumulative monitoring means from January 2019 through each method's first alarm with the 2016--2018 reference means. Quadratic contributions are signed terms $u_j(\\widehat V_m^{-1}u)_j$ divided by their sum; adjusted-range contributions are $u_j^2/\\widehat R_{m,j}^2$ divided by their sum.")
checkpoint_dates<-sort(unique(c("2019-12-06","2019-12-13","2019-12-31","2020-03-16","2020-03-31",amain$alarm_date,"2020-12-31")))
pp<-subset(paths,design=="continuous");rr<-subset(pp,method=="adjusted_range");qq<-subset(pp,method=="quadratic")
idx<-match(checkpoint_dates,rr$date)
tex_table("tables/spx_continuous_checkpoints.tex","tab:continuous-checkpoints","Selected dates under the single 2019--2020 monitoring horizon.","lrrr",
  "Date & Monitoring index & Adjusted-range ratio & Quadratic ratio",
  vapply(idx,function(i)sprintf("%s & %d & %.4f & %.4f",date_text(rr$date[i]),rr$k[i],rr$ratio[i],qq$ratio[i]),character(1)),
  "The monitoring index starts on January 2, 2019. Ratios use the two-year critical values throughout; the statistic is not restarted in January 2020.")
annual<-read_result("reference_annual")
tex_table("tables/spx_continuous_annual.tex","tab:continuous-annual","Annual summaries within the 2016--2018 reference sample.","lrcccc",
  "Year & $n$ & PC1 & PC2 & PC3 & PC4",
  vapply(unique(annual$year),function(y){x<-subset(annual,year==y);paste(c(y,x$n[1],sprintf("%.3f (%.3f)",x$mean/x$full_reference_sd,x$sd/x$full_reference_sd)),collapse=" & ")},character(1)),
  "Entries are the annual mean (standard deviation) of each original FPCA score, centered over 2016--2018 and divided by its full-reference standard deviation.")

png("figures/revised/spx_continuous_2019_2020.png",width=2100,height=1150,res=210)
layout(matrix(c(1,2),2,1),heights=c(1,.16));par(mar=c(4.8,4.5,1.5,1),mgp=c(2.6,.8,0),cex=1.1)
dd<-as.Date(rr$date);cols<-c("#1f5a99","#d47a17")
plot(dd,rr$ratio,type="l",col=cols[1],lwd=2,xlab="Date",ylab="Statistic / two-year 5% critical value",
     ylim=range(c(0,1,rr$ratio,qq$ratio)),xaxt="n")
lines(dd,qq$ratio,col=cols[2],lwd=2)
at<-as.Date(c("2019-01-02","2019-07-01","2020-01-02","2020-07-01","2020-12-31"))
axis.Date(1,at=at,format="%b %Y");abline(h=1,lty=2,col="black")
abline(v=as.numeric(as.Date("2020-01-01")),lty=3,col="grey65")
for(j in 1:2){k<-amain$alarm_index[j];points(dd[k],if(j==1)rr$ratio[k] else qq$ratio[k],pch=19,col=cols[j],cex=.9)}
par(mar=c(0,0,0,0));plot.new();legend("center",c("Adjusted-range SN","Quadratic SN","Alarm threshold"),
  col=c(cols,"black"),lty=c(1,1,2),lwd=2,horiz=TRUE,bty="n",text.col="black")
dev.off()
d<-readRDS("reproducibility/code/empirical/outputs/continuous_2019_2020/complete_analysis.rds")$continuous
png("figures/revised/spx_continuous_reference_scores.png",width=2000,height=1500,res=200)
par(mfrow=c(2,2),mar=c(4,4,2.5,1),cex=1)
zz<-sweep(sweep(d$ref,2,colMeans(d$ref),"-"),2,apply(d$ref,2,sd),"/")
for(j in 1:4) {plot(d$ref_dates,zz[,j],type="l",xlab="Date",ylab="Standardized reference score",main=paste("PC",j),
  col="#1f5a99",col.main="black");abline(h=0,lty=3,col="grey50")}
dev.off()
source("reproducibility/code/empirical/plot_spx_audit_figures.R")
ll<-read_result("loadings")
plot_spx_demixed_loadings(ll$u,as.matrix(ll[paste0("phi",1:4)]),as.matrix(ll[paste0("X",1:4)]),
 "figures/revised/spx_continuous_loadings.png")
