# Supplementary Figure S1. Sampling design; run from project root.
# Plot/subplot dimensions: Higuchi et al. (2012, 2015).
# C1 includes plants at least 20 cm tall.
# C1 occupies the central 10 m of the lower 0.5 m strip (nested sampling design).
if (.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','Portuguese_Brazil.utf8'))
library(grid)
out <- 'out/figures'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
s <- read.csv('data/sites.csv',fileEncoding='UTF-8')
d <- read.csv('data/multisitio_araucaria.csv',fileEncoding='UTF-8')
years <- lapply(split(d$ano,d$sitio),function(x) sort(unique(x)))
sites <- s[match(c('PB2','Alu','PARNAMUL','Epagri','Vargem','PAI','BJS'),s$code),]
sites$sampled_ha <- sites$parcelas*sites$area_parcela_m2/10000
sites$census_years <- vapply(sites$code,function(x) paste(years[[x]],collapse=', '),character(1))
stopifnot(sum(s$foco)==1, all(sites$sampled_ha>=.96 & sites$sampled_ha<=1),
 identical(years$PB1,c(2008L,2012L,2016L,2021L)))
W <- 10; H <- 7.4
ink <- '#253332'; grey <- '#637170'; line <- '#BAC6C2'
c1 <- '#467E79'; c2 <- '#A8C8C3'; c3 <- '#DDEBE7'; blue <- '#356C90'; red <- '#A74136'
txt <- function(x,y,label,size=10,face='plain',col=ink,just='left',rot=0) {
 grid.text(label,x=unit(x,'in'),y=unit(y,'in'),just=just,rot=rot,
  gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col,lineheight=1.15))
}
rect <- function(x,y,w,h,fill='white',col=line,lwd=.8) {
 grid.rect(x=unit(x,'in'),y=unit(y,'in'),width=unit(w,'in'),height=unit(h,'in'),
  just=c('left','bottom'),gp=gpar(fill=fill,col=col,lwd=lwd))
}
seg <- function(x,y,x2,y2,col=line,lwd=.8,arrow=NULL) {
 grid.lines(x=unit(c(x,x2),'in'),y=unit(c(y,y2),'in'),gp=gpar(col=col,lwd=lwd),arrow=arrow)
}
timeline <- function(x1,x2,y,yrs,col,labels=TRUE) {
 xx <- x1+(yrs-2008)/(2023-2008)*(x2-x1)
 seg(min(xx),y,max(xx),y,col,1.1)
 grid.points(x=unit(xx,'in'),y=unit(rep(y,length(xx)),'in'),pch=21,
  size=unit(.072,'in'),gp=gpar(fill=col,col=col))
 if(labels) for(i in seq_along(yrs)) txt(xx[i],y-.18,as.character(yrs[i]),size=8.7,just='center')
}
draw <- function(){
 grid.newpage()
 rect(0,0,W,H,fill='white',col=NA)
 seg(5.55,.35,5.55,7.06,col='#D9E0DC')
 txt(.3,7.04,'A  Focal population',14,'bold',red)
 txt(.3,6.72,'Pedras Brancas I',12,'bold')
 txt(.3,6.44,'50 permanent plots along six transects  |  1.00 ha sampled',10)
 txt(.3,6.10,'One 10 \u00d7 20 m plot',11,'bold')
 # Main plot is drawn to scale. The subplot strips share its lower edge.
 px <- .74; py <- 3.64; pw <- 4.30; ph <- 2.15
 rect(px,py,pw,ph,fill='#F5F7F6',col=ink,lwd=1.2)
 txt(px+pw/2,py+ph*.66,'Tree inventory',12,'bold',just='center')
 txt(px+pw/2,py+ph*.50,'DBH \u2265 5 cm  |  200 m\u00b2 per plot',10.5,just='center')
 txt(px+pw/2,py+.45,'Nested subplots: C1 within C2 within C3',9,col=grey,just='center')
 rect(px,py,pw,ph/10,c3,ink,.7)
 rect(px,py,pw,ph/20,c2,ink,.7)
 rect(px+pw*.25,py,pw/2,ph/20,c1,ink,.7)
 txt(px+pw*.50,py+ph/40,'C1',7.7,'bold','white',just='center')
 txt(px+pw*.125,py+ph/40,'C2',7.7,'bold',ink,just='center')
 txt(px+pw*.875,py+ph/40,'C2',7.7,'bold',ink,just='center')
 txt(px+pw*.5,py+ph*.075,'C3',7.7,'bold',ink,just='center')
 seg(px,py+ph+.11,px+pw,py+ph+.11,col=grey)
 txt(px+pw/2,py+ph+.23,'20 m',9,col=grey,just='center')
 seg(px-.14,py,px-.14,py+ph,col=grey)
 txt(px-.28,py+ph/2,'10 m',9,col=grey,just='center',rot=90)
 txt(.3,3.38,'Nested regeneration subplots (DBH < 5 cm)',10.5,'bold')
 txt(.3,3.08,'Class / height',9,'bold',grey)
 txt(2.44,3.08,'Dimensions / area',9,'bold',grey)
 txt(4.54,3.08,'Total\u00b9',9,'bold',grey,just='center')
 seg(.3,2.96,5.23,2.96)
 cl <- c('C1','C2','C3'); hs <- c('0.20 \u2264 h < 1 m','1 \u2264 h \u2264 3 m','h > 3 m')
 dims <- c('0.5 \u00d7 10 m / 5 m\u00b2','0.5 \u00d7 20 m / 10 m\u00b2','1 \u00d7 20 m / 20 m\u00b2')
 totals <- c('0.025 ha','0.050 ha','0.100 ha')
 for(i in 1:3){
  yy <- 2.77-(i-1)*.34
  rect(.3,yy-.07,.11,.14,c(c1,c2,c3)[i],col=grey,lwd=.5)
  txt(.48,yy,cl[i],10,'bold');txt(.91,yy,hs[i],9.5)
  txt(2.44,yy,dims[i],9.5);txt(4.54,yy,totals[i],9.5,just='center')
 }
 txt(.3,1.77,'\u00b9 50 subplots per class.  h, plant height.',8.7,col=grey)
 txt(.3,1.38,'Tree censuses',9.5,'bold')
 timeline(1.77,5.14,1.40,years$PB1,red)
 txt(.3,.83,'Regeneration censuses',9.5,'bold')
 timeline(1.77,5.14,.85,2011:2014,c1,labels=FALSE)
 txt(3.39,.85,'2011\u20132014 (annual)',9.3,col=c1)
 txt(.3,.32,expression(bold('Complete life-cycle model for ')*bolditalic('A. angustifolia')),10,col=red)
 # Regional sampling: explicit plot effort and actual inventory dates per site.
 txt(5.86,7.04,'B  Regional comparison',14,'bold',blue)
 txt(5.86,6.72,'Seven additional remnants',12,'bold')
 txt(5.86,6.44,expression(italic('A. angustifolia')*' trees with DBH \u2265 5 cm'),10.5)
 txt(5.86,6.13,'0.96\u20131.00 ha sampled per remnant',10)
 txt(5.86,5.75,'Remnant / censuses',9,'bold',grey)
 txt(8.38,5.75,'Plots \u00d7 m\u00b2',9,'bold',grey,just='center')
 txt(9.42,5.75,'Area (ha)',9,'bold',grey,just='center')
 seg(5.86,5.62,9.77,5.62)
 display <- c('Pedras Brancas II','Aluvial','Parque Natural Municipal','Amola Faca','Vargem','Painel','Bom Jardim da Serra')
 for(i in seq_len(nrow(sites))){
  yy <- 5.39-(i-1)*.52
  txt(5.86,yy,display[i],9.6,'bold')
  txt(5.86,yy-.18,sites$census_years[i],8.5,col=grey)
  txt(8.38,yy-.07,paste(sites$parcelas[i],'\u00d7',sites$area_parcela_m2[i]),9.5,just='center')
  txt(9.42,yy-.07,sprintf('%.2f',sites$sampled_ha[i]),9.5,just='center')
  seg(5.86,yy-.32,9.77,yy-.32,col='#E5EAE7',lwd=.6)
 }
 rect(5.86,.77,3.9,.90,fill='#EDF3F6',col=NA)
 txt(6.02,1.42,'Regional adult trajectories',10.5,'bold',blue)
 txt(6.02,1.12,'Focal population + seven comparison remnants',9)
 txt(6.02,.90,'Density and basal area',9)
 txt(5.86,.32,'DBH, diameter at breast height.',8.7,col=grey)
}
tiff(file.path(out,'Figure_S1.tiff'),width=W,height=H,units='in',res=600,compression='lzw',type='cairo');draw();dev.off()
cat('Figure S1 generated in ',out,'\n',sep='')
