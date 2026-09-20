suppressPackageStartupMessages(library(ggplot2))
dest <- file.path('out','figures')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
all <- as.data.frame(tables$comparison_summary)
all <- subset(all,scenario=='integrated_uniform')
metrics <- c('el_adult','el_stasis','el_advance','el_recruit')
rows <- match(metrics,all$metric)
stopifnot(!anyNA(rows))
elt <- all[rows,]
stopifnot(all(elt$lower<=elt$median),all(elt$median<=elt$upper))
elt$y <- 4:1
elt$label <- c('Adult\nsurvival','Stasis\n(C1\u2013C3)','Stage\nadvancement','Effective\nrecruitment')
lambda <- all[all$metric=='lambda',]
stopifnot(nrow(lambda)==1L)
annotation <- sprintf('\u03bb = %.3f [%.3f\u2013%.3f]',lambda$median,lambda$lower,lambda$upper)
p <- ggplot(elt,aes(x=median,y=y))+
  geom_col(orientation='y',fill='#4C72B0',width=.58,alpha=.90)+
  geom_segment(aes(x=lower,xend=upper,yend=y),linewidth=.6,colour='#303030')+
  geom_segment(aes(x=lower,xend=lower,y=y-.085,yend=y+.085),linewidth=.6,colour='#303030')+
  geom_segment(aes(x=upper,xend=upper,y=y-.085,yend=y+.085),linewidth=.6,colour='#303030')+
  geom_text(aes(y=y+.21,label=sprintf('%.2f',median)),hjust=-.2,size=3.4)+
  annotate('text',x=1.01,y=.49,label=annotation,hjust=1,size=3.4,colour='#454545')+
  scale_y_continuous(breaks=elt$y,labels=elt$label,limits=c(.38,4.5),expand=c(0,0))+
  scale_x_continuous(breaks=seq(0,1,.25),limits=c(0,1.04),expand=c(0,0))+
  labs(x='Elasticity of \u03bb',y=NULL)+
  theme_classic(base_size=11)+
  theme(panel.grid.major.y=element_line(colour='grey93',linewidth=.4),
        axis.ticks.y=element_line(colour='black',linewidth=.5),
        axis.line.y=element_line(colour='black',linewidth=.5),
        plot.margin=margin(7,8,7,7),axis.title.x=element_text(margin=margin(t=8)))
ggsave(file.path(dest,'Figure_3.tiff'),p,width=180,height=108,units='mm',dpi=600,device=grDevices::tiff,type='cairo',compression='lzw')
