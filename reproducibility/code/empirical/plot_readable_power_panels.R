# Replot saved scientific results; all five designs are retained in each figure.
p<-read.csv("reproducibility/data/paired_power_differences_5000rep.csv")
p<-subset(p,calibration=="size_adjusted" & change_shape=="abrupt")
designs<-c(gaussian_iid="DGP1: Gaussian IID",hong_var_diag="DGP2: Diagonal VAR",
           component_garch="DGP4: Component GARCH",component_sv="DGP5: Component SV",
           component_t5="DGP6: Student-t(5)")
for(direction in c("sparse_first","dense_first_two")){
  png(paste0("figures/revised/paired_power_abrupt_",direction,".png"),
      width=1800,height=1800,res=220,pointsize=15)
  layout(matrix(c(1:6,7,7),ncol=2,byrow=TRUE),heights=c(1,1,1,.22))
  par(mar=c(4.0,4.3,2.2,1),oma=c(0,0,0,0),cex=.95)
  for(id in names(designs)){
    d<-p[p$design==id & p$direction==direction,];d<-d[order(d$delta_standardized),]
    plot(d$delta_standardized,d$range_postchange_detection,type="b",pch=16,col="#0072B2",
         lwd=2,ylim=c(0,1),xlab="Standardized shift",ylab="Detection probability",
         main=designs[[id]],cex.main=1.0,xaxt="n",yaxt="n")
    axis(1,at=c(.05,.15,.30),labels=c("0.05","0.15","0.30"))
    axis(2,at=c(0,.5,1),labels=c("0.0","0.5","1.0"),las=1)
    abline(h=.05,lty=3,col="grey60")
    lines(d$delta_standardized,d$shao_postchange_detection,type="b",pch=17,col="#D55E00",lwd=2)
  }
  plot.new()
  par(mar=c(0,0,0,0))
  plot.new();legend("center",c("Adjusted-range SN","Quadratic SN","5% reference"),
                   col=c("#0072B2","#D55E00","grey60"),pch=c(16,17,NA),lty=c(1,1,3),
                   lwd=2,bty="n",cex=1.1,text.col="black",horiz=TRUE)
  dev.off()
}
