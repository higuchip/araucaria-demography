# Annual stage-structured model. Run from the repository root.
suppressPackageStartupMessages({library(data.table);library(popbio)})
root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
dest <- file.path(root, "out")
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
write_table <- function(x,name) assign(name,copy(x),envir=tables)
N <- 12000L
seed <- 20260912L
reg<-fread(file.path(root,"data","focal_regeneracao.csv"),encoding="UTF-8",na.strings=c("NA",""))
ad<-fread(file.path(root,"data","focal_arboreo.csv"),encoding="UTF-8")
reg<-reg[grepl("^Araucaria angustifolia",especie)]
ad<-ad[especie=="Araucaria angustifolia"]
reg[,id:=paste(transecto,parcela,subparcela,individuo,sep="_")]
ad[,id:=paste(transecto,parcela,arvore,sep="_")]
stopifnot(!anyDuplicated(reg[,.(id,ano)]),!anyDuplicated(ad[,.(id,ano)]))
setorder(reg,id,ano);reg[,next_class:=shift(classe,-1),by=id];reg[,next_year:=shift(ano,-1),by=id]
tr<-reg[!is.na(next_year)&next_year==ano+1L]
cnt<-function(from,to)nrow(tr[classe==from&next_class==to])
n1<-c(C1=cnt("1","1"),C2=cnt("1","2"),M=cnt("1","M"))
n2<-c(C1=cnt("2","1"),C2=cnt("2","2"),C3=cnt("2","3"),M=cnt("2","M"))
n3<-c(C2=cnt("3","2"),C3=cnt("3","3"),M=cnt("3","M"))
k2<-cnt("2","AD");k3<-cnt("3","AD");E2<-sum(n2)+k2;E3<-sum(n3)+k3
stopifnot(k2==1L,k3==0L,E2==11L,E3==26L)
areas<-c("1"=.025,"2"=.05,"3"=.1,"A"=1) # hectares, entire sampled network
dens<-reg[classe%in%c("1","2","3"),.(n=.N),by=.(ano,classe)]
dens[,density:=n/areas[classe]]
D2all<-mean(dens[classe=="2",density]);D3all<-mean(dens[classe=="3",density])
nyreg<-diff(range(reg$ano));D2start<-E2/(nyreg*areas['2']);D3start<-E3/(nyreg*areas['3'])
adult_density<-ad[status%in%c("viva","recrutamento"),.N,by=ano][order(ano)]
DA<-mean(adult_density$N)/areas['A']
K<-ad[status=="recrutamento",.N];Tadult<-diff(range(ad$ano));Krest<-K-k2-k3
stopifnot(K==11L,Tadult==13L,Krest==10L)
recN<-tr[is.na(classe)&next_class=="1",.N]
# Fail on previously unhandled live-origin destinations instead of dropping them.
stopifnot(sum(n1)==nrow(tr[classe=="1"]),E2==nrow(tr[classe=="2"]),E3==nrow(tr[classe=="3"]))
write_table(tr[,.N,by=.(from=classe,to=next_class)],"observed_fates")
write_table(dens,"regeneration_densities")
write_table(data.table(component=c("C1","C2","C3","A"),area_ha=as.numeric(areas),
  start= c(rep(min(reg$ano),3),min(ad$ano)),end=c(rep(max(reg$ano),3),max(ad$ano))),"sampling")
write_table(data.table(quantity=c("K_total","K_known_C2","K_known_C3","K_complement","T_adult","E2_individual_years","E3_individual_years","D2_all","D3_all","D2_interval_start","D3_interval_start","C1_recruits","mean_A_density"),
 value=c(K,k2,k3,Krest,Tadult,E2,E3,D2all,D3all,D2start,D3start,recN,DA)),"integration_inputs")
aw<-dcast(ad,id~ano,value.var="status")
yrs<-sort(unique(ad$ano));iv<-lapply(seq_len(length(yrs)-1L),function(i)c(yrs[i],yrs[i+1]))
risk<-function(include_recruits)rbindlist(lapply(iv,function(v){
 start<-aw[[as.character(v[1])]]%in%if(include_recruits)c("viva","recrutamento")else"viva"
 end<-aw[[as.character(v[2])]][start]
 stopifnot(!anyNA(end),all(end%in%c("viva","morta","recrutamento")))
 data.table(start=v[1],end=v[2],years=diff(v),alive=sum(end%in%c("viva","recrutamento")),dead=sum(end=="morta"))
}))
rn<-risk(TRUE);write_table(rn,"adult_risk_sets")
set.seed(seed)
rdirn<-function(n,alpha){m<-matrix(rgamma(n*length(alpha),rep(alpha,each=n)),nrow=n);m/rowSums(m)}
p1<-rdirn(N,n1+.5);p2<-rdirn(N,n2+.5);p3<-rdirn(N,n3+.5)
uS<-matrix(runif(N*nrow(rn)),N)
SN<-sapply(1:nrow(rn),function(j)qbeta(uS[,j],rn$alive[j]+.5,rn$dead[j]+.5)^(1/rn$years[j]))
sNew<-rowMeans(SN)
birth<-rgamma(N,recN+.5,rate=nyreg)/areas['1']/DA
# Preserve the random-number sequence of the reported posterior sample.
invisible(rgamma(N,K+.5,rate=Tadult*areas['A']))
uSample<-runif(N)

# Approximate joint observation model:
# q2=rho*I/D2, q3=(1-rho)*I/D3; I is effective entries ha^-1 yr^-1.
# Known sampled exits use Binomial(Ej,qj), non-A fates use conditional multinomials.
# The remaining TEN entries have Poisson mean I*Tadult*area_A-E2*q2-E3*q3.
# This uses the confirmed overlapping event once; the ten have unknown stage origins.
# Priors: density(I) proportional to I^-1/2; rho ~ Beta(a,b), independent.
# Proposal I~Gamma(K+.5,Tadult*area_A), rho~Beta(a+k2,b+k3).
integration_draws<-function(a,b,D2,D3,tag,propN=200000L,extra_seed=0L){
 set.seed(seed+100L+extra_seed)
 Ip<-rgamma(propN,K+.5,rate=Tadult*areas['A']);rhop<-rbeta(propN,a+k2,b+k3)
 q2p<-rhop*Ip/D2;q3p<-(1-rhop)*Ip/D3
 H<-Tadult*areas['A']-E2*rhop/D2-E3*(1-rhop)/D3
 ok<-q2p>0&q2p<1&q3p>0&q3p<1&H>0
 logw<-rep(-Inf,propN)
 logw[ok]<-Krest*log(H[ok]/(Tadult*areas['A']))+Ip[ok]*(Tadult*areas['A']-H[ok])+
   (E2-k2)*log1p(-q2p[ok])+(E3-k3)*log1p(-q3p[ok])
 w<-exp(logw-max(logw));w<-w/sum(w);ESS<-1/sum(w*w)
 # Independent algebra check against full target / proposal log density.
 ix<-which(ok)[seq_len(100L)]
 lt<-dbinom(k2,E2,q2p[ix],log=TRUE)+dbinom(k3,E3,q3p[ix],log=TRUE)+
     dpois(Krest,Ip[ix]*H[ix],log=TRUE)-.5*log(Ip[ix])+dbeta(rhop[ix],a,b,log=TRUE)
 lp<-dgamma(Ip[ix],K+.5,rate=Tadult*areas['A'],log=TRUE)+dbeta(rhop[ix],a+k2,b+k3,log=TRUE)
 residual<-(lt-lp)-logw[ix]
 stopifnot(diff(range(residual))<1e-8,ESS>.8*propN,max(w)<.001)
 pick<-findInterval(uSample,c(0,cumsum(w)),all.inside=TRUE)
 means<-c(I=sum(w*Ip),rho=sum(w*rhop))
 batches<-sapply(split(seq_len(propN),rep(1:20,length.out=propN)),function(ii){
    wi<-w[ii]/sum(w[ii]);c(sum(wi*Ip[ii]),sum(wi*rhop[ii]))})
 diagnostics<-data.table(scenario=tag,proposals=propN,ESS=ESS,ESS_fraction=ESS/propN,
  rejected=sum(!ok),weight_max=max(w),I_mean=means['I'],rho_mean=means['rho'],
  I_MCSE=sd(batches[1,])/sqrt(20),rho_MCSE=sd(batches[2,])/sqrt(20),density_ratio_error=diff(range(residual)))
 list(I=Ip[pick],rho=rhop[pick],q2=q2p[pick],q3=q3p[pick],D2=D2,D3=D3,diagnostics=diagnostics)
}
fits<-list(integrated_uniform=integration_draws(1,1,D2start,D3start,"integrated_uniform"))
write_table(rbindlist(lapply(fits,`[[`,"diagnostics")),"importance_diagnostics")

build_matrix<-function(i,mode,fit=NULL){
 A<-matrix(0,4,4,dimnames=list(c("C1","C2","C3","A"),c("C1","C2","C3","A")))
 A[1:2,1]<-p1[i,1:2];A[4,4]<-sNew[i];A[1,4]<-birth[i]
 q2<-fit$q2[i];q3<-fit$q3[i]
 A[1:3,2]<-(1-q2)*p2[i,1:3];A[4,2]<-q2
 A[c(2,3),3]<-(1-q3)*p3[i,1:2];A[4,3]<-q3
 A
}
eigen_lambda<-function(A)max(Re(eigen(A,only.values=TRUE)$values))
matrix_metrics<-function(A){
 ea<-eigen.analysis(A);el<-ea$elasticities;F<-A*0;F[1,4]<-A[1,4];U<-A-F
 stopifnot(all(A>=0),all(colSums(U)<=1+1e-10),abs(sum(el)-1)<1e-8)
 NN<-solve(diag(4)-U);R0<-eigen_lambda(F%*%NN)
 NP<-solve(diag(3)-U[1:3,1:3]);h<-as.numeric(U[4,1:3]%*%NP[,1])
 # At lambda=1: required annual adult death probability = effective recruitment * P(reach A).
 qstar<-A[1,4]*h
 c(lambda=ea$lambda1,R0=R0,T_R0=log(R0)/log(ea$lambda1),T_age_diff=1/el[1,4],
  adult_residence=1/-log(A[4,4]),adult_expected_years=1/(1-A[4,4]),coast10=log(.1)/log(A[4,4]),coast50=log(.5)/log(A[4,4]),C1_residence=NN[1,1],
  reach_A_probability=h,conditional_years_to_A=as.numeric(U[4,1:3]%*%(NP%*%NP)[,1])/h,
  el_adult=el[4,4],el_stasis=sum(diag(el)[1:3]),el_advance=el[2,1]+el[3,2]+el[4,2]+el[4,3],
  el_recruit=el[1,4],el_retro=el[1,2]+el[2,3],damping=ea$damping.ratio,
  q2A=A[4,2],q3A=A[4,3],g12=A[2,1],g23=A[3,2],sA=A[4,4],b=A[1,4],mortality_threshold=qstar)
}
scenarios<-names(fits)
all_results<-list();matrix_draws<-list()
for(s in scenarios){
 cat("Calculating",s,"\n");flush.console()
 fit<-fits[[s]];mm<-vector("list",N);Amat<-array(NA_real_,c(4,4,N),dimnames=list(c("C1","C2","C3","A"),c("C1","C2","C3","A"),NULL))
 for(i in seq_len(N)){
  A<-build_matrix(i,s,fit);Amat[,,i]<-A;mm[[i]]<-matrix_metrics(A)
  if(i<=12L){
    F<-A*0;F[1,4]<-A[1,4];U<-A-F
    stopifnot(abs(Rage::gen_time(U,F,method="R0")-mm[[i]]['T_R0'])<1e-7)
    stopifnot(abs(Rage::gen_time(U,F,method="age_diff")-mm[[i]]['T_age_diff'])<1e-7)
    stopifnot(abs(Rage::net_repro_rate(U,F)-mm[[i]]['R0'])<1e-7)
    if(mm[[i]]['mortality_threshold']<1){
      B<-A;B[4,4]<-1-mm[[i]]['mortality_threshold'];stopifnot(abs(eigen_lambda(B)-1)<1e-8)
    }
    B<-A;B[4,4]<-A[4,4]*1.000001
    fd<-(eigen_lambda(B)-mm[[i]]['lambda'])/mm[[i]]['lambda']/1e-6
    stopifnot(abs(fd-mm[[i]]['el_adult'])<1e-5)
  }
 }
 tab<-as.data.table(do.call(rbind,mm));tab[,`:=`(draw=1:N,scenario=s)]
 D2<-if(is.null(fit))D2all else fit$D2;D3<-if(is.null(fit))D3all else fit$D3
 tab[,ingrowth_intensity:=q2A*D2+q3A*D3]
 tab[,expected_total_entries:=Tadult*areas['A']*ingrowth_intensity]
 tab[,expected_sample_C2:=E2*q2A];tab[,expected_sample_C3:=E3*q3A]
 if(!is.null(fit)){
  stopifnot(max(abs(tab$ingrowth_intensity-fit$I))<1e-12)
  tab[,rho:=fit$rho]
 }else tab[,rho:=ifelse(ingrowth_intensity>0,q2A*D2/ingrowth_intensity,NA_real_)]
 # Expected totals decompose into disjoint observed exposures; no duplicated event.
 tab[,expected_complement:=expected_total_entries-expected_sample_C2-expected_sample_C3]
 stopifnot(all(tab$expected_complement>=0),max(abs(tab$expected_total_entries-
   tab$expected_complement-tab$expected_sample_C2-tab$expected_sample_C3))<1e-10)
 all_results[[s]]<-tab;matrix_draws[[s]]<-Amat
}
draws<-rbindlist(all_results,use.names=TRUE,fill=TRUE)
tables$comparison_draws <- draws
summ<-melt(draws,id.vars=c("scenario","draw"),variable.name="metric",value.name="value")
summ<-summ[,.(median=median(value),lower=quantile(value,.025),upper=quantile(value,.975),n=.N),by=.(scenario,metric)]
write_table(summ,"comparison_summary")
probs<-draws[,.(P_lambda_ge1=mean(lambda>=1),P_R0_gt1=mean(R0>1),
 P_adult_largest_group=mean(el_adult>pmax(el_stasis,el_advance,el_recruit,el_retro)),
 P_threshold_gt1=mean(mortality_threshold>1),
 P_C2_at_least_one=mean(1-(1-q2A)^E2),P_C3_zero=mean((1-q3A)^E3)),by=scenario]
write_table(probs,"comparison_probabilities")
# Single preserved set of matrices for all downstream metrics and figures.
tables$annual_model <- list(matrices=matrix_draws,draws=draws,adult_survival_intervals=SN,
 risk_sets=rn,fits=fits,seed=seed,N=N)
source(file.path(root,"R","02_scenarios_transients.R"),local=TRUE)
run_derived_analyses(matrix_draws,all_results,SN,rn,dens,adult_density,DA,root,dest,tables)

# Posterior predictive checks: same observation process, known exits plus complement.
set.seed(seed+300L)
ppc<-rbindlist(lapply(scenarios,function(s){
 x<-all_results[[s]];y2<-rbinom(N,E2,x$q2A);y3<-rbinom(N,E3,x$q3A)
 yr<-rpois(N,x$expected_complement);yt<-y2+y3+yr
 data.table(scenario=s,measure=c("C2_to_A","C3_to_A","complement","total"),
 observed=c(k2,k3,Krest,K),
 predicted_median=c(median(y2),median(y3),median(yr),median(yt)),
 predicted_lower=c(quantile(y2,.025),quantile(y3,.025),quantile(yr,.025),quantile(yt,.025)),
 predicted_upper=c(quantile(y2,.975),quantile(y3,.975),quantile(yr,.975),quantile(yt,.975)),
 P_rep_ge_observed=c(mean(y2>=k2),mean(y3>=k3),mean(yr>=Krest),mean(yt>=K)),
 P_rep_le_observed=c(mean(y2<=k2),mean(y3<=k3),mean(yr<=Krest),mean(yt<=K)))
}))
write_table(ppc,"posterior_predictive_checks")
metrics_key<-c("lambda","el_adult","R0","T_R0","q2A","q3A","rho","ingrowth_intensity")
write_table(summ[metric%in%metrics_key],"key_results")
cat("Annual model completed.\n")
print(summ[metric %in% c("lambda","R0","T_R0","el_adult")])
print(probs)
