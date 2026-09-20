suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
dest <- file.path('out','figures')
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
lab <- c('5-15','15-25','25-35','35-45','45+')
aa <- copy(tables$ad_ara)
aa <- aa[status %in% c('viva','recrutamento') & ano %in% c(2008L,2021L)]
aa[, cls := cut(dap_cm, breaks=c(5,15,25,35,45,Inf), labels=lab, right=FALSE)]
cnt <- aa[, .(n=.N), by=.(cls=factor(cls, levels=lab), ano=factor(ano))]
d <- merge(CJ(cls=factor(lab, levels=lab), ano=factor(c('2008','2021'))),
           cnt, by=c('cls','ano'), all.x=TRUE)
d[is.na(n), n := 0L]
sampled_area_ha <- 1
d[, density := n / sampled_area_ha]
stopifnot(identical(as.integer(d[ano=='2008', n]), c(66L,25L,17L,15L,7L)),
          identical(as.integer(d[ano=='2021', n]), c(56L,24L,14L,16L,21L)))
p <- ggplot(d, aes(cls, density, fill=ano)) +
  geom_col(position=position_dodge(.7), width=.62) +
  geom_text(aes(label=n), position=position_dodge(.7), vjust=-.4, size=3.2) +
  scale_fill_manual(values=c('2008'='#B5C7D3','2021'='#2E6E8E')) +
  labs(x='DBH class (cm)',
       y=expression('Adult araucaria density (individuals ha'^-1*')'), fill=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position='inside', legend.position.inside=c(.85,.85),
        panel.grid.major.y=element_line(colour='grey92'))
ggsave(file.path(dest,'Figure_5.tiff'),p,width=6.4,height=4.2,units='in',dpi=600,device=grDevices::tiff,type='cairo',compression='lzw')
