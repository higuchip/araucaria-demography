# All analyses consume the exact matrices calculated by 01_annual_model.R.
run_derived_analyses<-function(matrices,results,SN,risk,dens,adult_density,DA,root,dest,tables){
 write_table<-function(x,name)assign(name,data.table::copy(x),envir=tables)
 lam<-function(A)max(Re(eigen(A,only.values=TRUE)$values))
 qs<-function(z)c(median=median(z),lower=unname(quantile(z,.025)),upper=unname(quantile(z,.975)))
 N<-dim(matrices[[1]])[3]
 primary<-"integrated_uniform"
 
 # Four-state matrices and complete survival fates, including death.
 cells<-rbindlist(lapply(names(matrices),function(s)rbindlist(lapply(1:4,function(j)rbindlist(lapply(1:4,function(i){
  z<-matrices[[s]][i,j,];v<-qs(z)
  data.table(model=s,from=dimnames(matrices[[s]])[[2]][j],to=dimnames(matrices[[s]])[[1]][i],
   median=v[1],lower=v[2],upper=v[3])
 }))))))
 write_table(cells,"matrix_elements")
 deaths<-rbindlist(lapply(names(matrices),function(s)rbindlist(lapply(1:4,function(j){
  z<-1-colSums(matrices[[s]][,j,]);if(j==4)z<-1-matrices[[s]][4,4,]
  stopifnot(min(z)>-1e-10);v<-qs(z)
  data.table(model=s,from=dimnames(matrices[[s]])[[2]][j],median=v[1],lower=v[2],upper=v[3])
 }))))
 write_table(deaths,"annual_death_probabilities")

 # Survival annualization is checked at the observed duration of each interval.
 stopifnot(nrow(SN)==N,ncol(SN)==nrow(risk))
 obs<-copy(risk)
 obs[,`:=`(S_observed=alive/(alive+dead),s_annual=(alive/(alive+dead))^(1/years))]
 obs[,`:=`(q_annual=1-s_annual,mu_equivalent=-log(s_annual))]
 stopifnot(max(abs(obs$s_annual^obs$years-obs$S_observed))<1e-12)
 obs[,interval:=paste(start,end,sep="-")]
 write_table(obs,"observed_mortality_intervals")
 mort<-rbindlist(lapply(seq_len(nrow(obs)),function(j){
  q<-1-SN[,j];v<-qs(q)
  data.table(interval=obs$interval[j],q_median=v[1],q_lower=v[2],q_upper=v[3],
   q_observed=obs$q_annual[j],years=obs$years[j])
 }))
 write_table(mort,"mortality_interval_posterior")
 
 # Each interval is held constant indefinitely; these are not finite mortality pulses.
 interval_draws<-rbindlist(lapply(names(matrices),function(s){
  cat("Interval scenarios:",s,"\n");flush.console()
  rbindlist(lapply(seq_len(nrow(obs)),function(j){
   z<-vapply(seq_len(N),function(i){A<-matrices[[s]][,,i];A[4,4]<-SN[i,j];lam(A)},numeric(1))
   data.table(model=s,interval=obs$interval[j],draw=seq_len(N),lambda=z,
    annual_death=1-SN[,j],delta_lambda_vs_reference=z-results[[s]]$lambda,
    q_above_threshold=(1-SN[,j])>=results[[s]]$mortality_threshold)
  }))
 }))
 write_table(interval_draws[,.(lambda_median=median(lambda),lambda_lower=quantile(lambda,.025),
 lambda_upper=quantile(lambda,.975),P_lambda_ge1=mean(lambda>=1),
 delta_median=median(delta_lambda_vs_reference),delta_lower=quantile(delta_lambda_vs_reference,.025),
 delta_upper=quantile(delta_lambda_vs_reference,.975),P_mortality_ge_threshold=mean(q_above_threshold)),
 by=.(model,interval)],"interval_scenarios")
 stopifnot(all(interval_draws$q_above_threshold==(interval_draws$lambda<=1)))
 tables$interval_scenario_draws <- interval_draws

 # Shared-nuisance contrast: higher mortality interval minus 2012-2016.
 contrasts<-rbindlist(lapply(names(matrices),function(s){
  a<-interval_draws[model==s & interval=="2012-2016"][order(draw)]
  b<-interval_draws[model==s & interval=="2016-2021"][order(draw)]
  z<-b$lambda-a$lambda;v<-qs(z)
  data.table(model=s,contrast="2016-2021 minus 2012-2016",median=v[1],lower=v[2],upper=v[3],P_negative=mean(z<0))
 }))
 write_table(contrasts,"mortality_scenario_contrasts")

 # Primary annual curve; q is a probability, not an instantaneous hazard.
 grid<-seq(0,.1,.001)
 L<-matrix(NA_real_,N,length(grid))
 for(k in seq_along(grid)){
  L[,k]<-vapply(seq_len(N),function(i){A<-matrices[[primary]][,,i];A[4,4]<-1-grid[k];lam(A)},numeric(1))
  if(k%%20L==1L){cat("Mortality curve",k,"/",length(grid),"\n");flush.console()}
 }
 stopifnot(all(L[,-1,drop=FALSE]-L[,-ncol(L),drop=FALSE]<=1e-10))
 curve<-data.table(q_annual=grid,median=apply(L,2,median),lower=apply(L,2,quantile,.025),
 upper=apply(L,2,quantile,.975),P_lambda_ge1=colMeans(L>=1))
 write_table(curve,"annual_mortality_curve")
 threshold<-rbindlist(lapply(names(results),function(s){
  q<-results[[s]]$mortality_threshold;v<-qs(q);ok<-q<1
  muv<-if(all(ok))qs(-log1p(-q))else rep(NA_real_,3)
  data.table(model=s,q_median=v[1],q_lower=v[2],q_upper=v[3],P_no_root_in_0_1=mean(!ok),
   mu_equivalent_median=muv[1],mu_equivalent_lower=muv[2],mu_equivalent_upper=muv[3])
 }))
 write_table(threshold,"annual_mortality_threshold")
 crossing<-approx(curve$median,curve$q_annual,xout=1)$y
 write_table(data.table(quantity=c("median_of_draw_specific_thresholds","crossing_of_pointwise_median_curve"),
 value=c(median(results[[primary]]$mortality_threshold),crossing)),"threshold_curve_crossing")
 # Initial stage densities from the concurrent 2012 census.
 year<-2012L
 n2012<-c(vapply(c("1","2","3"),function(cl)dens[ano==year & classe==cl,density],numeric(1)),
 adult_density[ano==year,N]/1)
 names(n2012)<-c("C1","C2","C3","A")
 n0s<-list(concurrent_2012=n2012)
 write_table(rbindlist(lapply(names(n0s),function(k)data.table(vector=k,stage=names(n0s[[k]]),density_ha=n0s[[k]]))),"transient_initial_vectors")
 
 project<-function(M,lambda,n0){
  x<-matrix(n0/sum(n0),4,N);best<-best50<-best20<-rep(1,N);peak<-rep(0L,N)
  old_best<-rep(-Inf,N);minimum<-rep(1,N);values<-list();checks<-integer(0)
  for(t in 1:400){
   nx<-matrix(0,4,N)
   for(j in 1:4)for(i in 1:4)nx[i,]<-nx[i,]+M[i,j,]*x[j,]
   x<-sweep(nx,2,lambda,"/");total<-colSums(x)
   stopifnot(all(is.finite(total)),all(total>=0))
   improved<-total>best+1e-12;peak[improved]<-t;best<-pmax(best,total)
   old_best<-pmax(old_best,total);minimum<-pmin(minimum,total)
   if(t<=20)best20<-pmax(best20,total)
   if(t<=50)best50<-pmax(best50,total)
   if(t%in%c(1,5,10,20,50,100,200,400)){
    values[[as.character(t)]]<-total
    if(t==5){
     for(k in 1:min(12,N)){
      checkx<-n0/sum(n0)
      for(tt in 1:5)checkx<-(M[,,k]/lambda[k])%*%checkx
      stopifnot(max(abs(as.numeric(checkx)-x[,k]))<1e-10)
     }
    }
   }
  }
  stopifnot(all(best>=1),all(best50<=best+1e-10),all(best20<=best50+1e-10))
  list(draws=data.table(draw=seq_len(N),maxamp_0_400=best,maxamp_0_50=best50,maxamp_0_20=best20,
   peak_year_0_400=peak,maxamp_1_400_legacy_definition=old_best,minamp_0_400=minimum,relative_total_400=total),
   trajectory=rbindlist(lapply(names(values),function(t){v<-qs(values[[t]]);data.table(year=as.integer(t),median=v[1],lower=v[2],upper=v[3])})))
 }
 transients<-list();trajectories<-list()
 for(s in names(matrices))for(k in names(n0s)){
  cat("Transient projection:",s,k,"\n");flush.console()
  tr<-project(matrices[[s]],results[[s]]$lambda,n0s[[k]])
  tr$draws[,`:=`(model=s,initial_vector=k)];tr$trajectory[,`:=`(model=s,initial_vector=k)]
  key<-paste(s,k,sep="__");transients[[key]]<-tr$draws;trajectories[[key]]<-tr$trajectory
 }
 td<-rbindlist(transients)
 td[,peak_year_0_400:=as.numeric(peak_year_0_400)]
 tables$transient_draws <- td
 tall<-melt(td,id.vars=c("model","initial_vector","draw"),variable.name="metric")
 write_table(tall[,.(median=median(value),lower=quantile(value,.025),upper=quantile(value,.975)),
 by=.(model,initial_vector,metric)],"transient_summary")
 write_table(td[,.(P_peak_at_year_400=mean(peak_year_0_400==400),P_peak_at_year_0=mean(peak_year_0_400==0),
 P_maxamp_exceeds_1_1=mean(maxamp_0_400>1.1)),by=.(model,initial_vector)],"transient_diagnostics")
 write_table(rbindlist(trajectories),"transient_trajectory_summary")

 write_table(data.table(check=c("annual_survival_reconstructs_observed_intervals","mortality_curve_monotone", "interval_lambda_threshold_agree",
 "transient_batched_vs_individual_projection","transient_includes_initial_time","shared_matrices_for_all_analyses"),passed=TRUE),"derived_verification")
 invisible(TRUE)
}
