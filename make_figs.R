source("engine2.R")
set.seed(11L)
n_r0<-250L; n_ec<-500L; Delta<-4; piB<-0.30
dat<-gen_data(1L,2L,n_r0,n_ec,"A",Delta,piB,FALSE)
Yc<-dat$Yc;Xc<-dat$Xc;Ye<-dat$Ye;Xe<-dat$Xe;B<-dat$B;cc<-C2
Zr<-as.numeric(Yc%*%cc);Ze<-as.numeric(Ye%*%cc)
mf<-fit_mean(Yc,Xc,"correct");covR<-cov_stratum(Yc-mf(Xc),Xc)
scr<-score_vec(Yc-mf(Xc),Xc,covR,"c2");sce<-score_vec(Ye-mf(Xe),Xe,covR,"c2")
q<-as.numeric(quantile(abs(scr),0.85))
gray<-"#888780";blue<-"#378ADD";amber<-"#EF9F27";green<-"#97C459"
tr<-function(c,a=0.32)adjustcolor(c,alpha.f=a)
fc<-mean(B==0);fk<-mean(B==1)
dc<-function(v,sc=1){d<-density(v,n=220,adjust=1.1);d$y<-d$y*sc;d}
root<-""   # write to current working directory (run from repo root)
lg<-function() legend("topright",bty="n",cex=0.85,legend=c("RCT control","compatible EC  (borrow)","non-compatible EC  (drop)"),fill=c(tr(gray),tr(blue),tr(amber)),border=c(gray,blue,amber))
# FIG 1 : outcome Y
pdf(paste0(root,"fig_outcome.pdf"),width=6.6,height=3.0,pointsize=10)
par(mar=c(3.2,0.8,0.6,0.6),mgp=c(1.9,0.6,0))
d1<-dc(Zr);d2<-dc(Ze[B==0],fc);d3<-dc(Ze[B==1],fk)
plot(NA,xlim=range(d1$x,d2$x,d3$x),ylim=c(0,max(d1$y)),xlab="outcome Y  (visit-average control level)",ylab="",yaxt="n",bty="n")
polygon(d3,col=tr(amber),border=amber,lwd=1.6);polygon(d2,col=tr(blue),border=blue,lwd=1.6);polygon(d1,col=tr(gray),border=gray,lwd=1.6);lg()
dev.off()
# FIG 2 : nonconformity score
pdf(paste0(root,"fig_score.pdf"),width=6.6,height=3.0,pointsize=10)
par(mar=c(3.2,0.8,0.6,0.6),mgp=c(1.9,0.6,0))
d1<-dc(scr);d2<-dc(sce[B==0],fc);d3<-dc(sce[B==1],fk);ym<-max(d1$y)*1.1
plot(NA,xlim=range(d1$x,d2$x,d3$x),ylim=c(0,ym),xlab="nonconformity score  (standardized residual)",ylab="",yaxt="n",bty="n")
rect(-q,0,q,ym,col=tr(green,0.16),border=NA);abline(v=c(-q,q),col="#639922",lty=2)
polygon(d3,col=tr(amber),border=amber,lwd=1.6);polygon(d2,col=tr(blue),border=blue,lwd=1.6);polygon(d1,col=tr(gray),border=gray,lwd=1.6)
text(0,ym*0.97,"sym-ada keeps this band",col="#3B6D11",cex=0.82)
text(mean(sce[B==1]),ym*0.42,"drop",col="#A32D2D",cex=0.9);lg()
dev.off()
cat("done q=",round(q,2)," Ymeans",round(mean(Zr),2),round(mean(Ze[B==0]),2),round(mean(Ze[B==1]),2)," Smeans",round(mean(scr),2),round(mean(sce[B==0]),2),round(mean(sce[B==1]),2),"\n")
