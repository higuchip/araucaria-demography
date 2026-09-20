suppressPackageStartupMessages(library(ggplot2))
src <- file.path('out')
dest <- file.path('out','figures')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
cur <- as.data.frame(tables$annual_mortality_curve)
obs <- as.data.frame(tables$observed_mortality_intervals)
thresholds <- as.data.frame(tables$annual_mortality_threshold)
threshold <- subset(thresholds,model=='integrated_uniform')
stopifnot(nrow(threshold)==1L,nrow(obs)==3L,all(diff(cur$median)<=0),
          all(cur$lower<=cur$median),all(cur$median<=cur$upper))
root <- 100*threshold$q_median
obs$percent <- 100*obs$q_annual
obs$label <- c('2008\u20132012 (0%)','2012\u20132016 (0.38%)','2016\u20132021 (1.20%)')
obs$label_y <- c(1.066,1.056,1.046)
obs$colour <- c('#606060','#0F6E56','#98700A')
plot <- ggplot(cur,aes(x=100*q_annual,y=median))+
  annotate('rect',xmin=0,xmax=10,ymin=-Inf,ymax=1,fill='#C44E52',alpha=.05)+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill='#DD8452',alpha=.22)+
  geom_hline(yintercept=1,linetype='dashed',colour='grey40')+
  geom_line(linewidth=1.1,colour='#B2480F')+
  geom_vline(data=obs,aes(xintercept=percent,colour=colour),linetype='dotted',linewidth=.65,show.legend=FALSE)+
  scale_colour_identity()+
  geom_vline(xintercept=root,colour='#C44E52',linewidth=.7)+
  geom_label(data=obs,aes(x=percent+.12,y=label_y,label=label,colour=colour),
             inherit.aes=FALSE,hjust=0,size=2.9,fill='white',linewidth=0,
             label.padding=grid::unit(.10,'lines'),show.legend=FALSE)+
  annotate('text',x=root+.15,y=.948,hjust=0,size=3.1,colour='#8C2F31',
           label=sprintf('Median threshold: %.2f%%',root))+
  annotate('text',x=9.65,y=.948,hjust=1,size=3.2,colour='#8C2F31',label='Decline')+
  scale_x_continuous(breaks=seq(0,10,2.5),limits=c(0,10),expand=expansion(add=c(.14,0)))+
  scale_y_continuous(breaks=c(.94,.98,1,1.02,1.06),expand=expansion(mult=c(.025,.04)))+
  labs(x='Annual adult mortality (%)',y=expression('Population growth ('*lambda*')'))+
  theme_classic(base_size=11)+
  theme(panel.grid.major.y=element_line(colour='grey92'),plot.margin=margin(8,10,7,7),
        axis.line=element_line(colour='black',linewidth=.5))
ggsave(file.path(dest,'Figure_4.tiff'),plot,width=180,height=120,units='mm',dpi=600,device=grDevices::tiff,type='cairo',compression='lzw')
