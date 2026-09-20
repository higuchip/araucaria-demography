# Figure 1: municipalities, study remnants, and original Araucaria Forest.
# Run from repository root. Third-party layers are downloaded to the R user cache.
.libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths()))
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE','Portuguese_Brazil.utf8'))
suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(ggspatial); library(cowplot)
})
sf_use_s2(FALSE)
options(timeout=max(600,getOption('timeout')))
root <- normalizePath('.', winslash='/')
out <- file.path(root, 'out/figures')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
map_cache <- getOption('araucaria.map_cache', tools::R_user_dir('araucaria-demography','cache'))
cache <- file.path(map_cache,'ibge_vegetation_2004')
dir.create(cache, recursive=TRUE, showWarnings=FALSE)
url <- 'https://geoftp.ibge.gov.br/informacoes_ambientais/vegetacao/vetores/brasil_5000_mil/'
zip <- file.path(cache, 'Vegetacao_5000mil.zip')
if (!file.exists(zip)) download.file(paste0(url, 'Vegetacao_5000mil.zip'), zip, mode='wb')
if (!file.exists(file.path(cache, 'Vegetacao_5000.shp'))) unzip(zip, exdir=cache)
v <- st_read(file.path(cache, 'Vegetacao_5000.shp'), options='ENCODING=UTF-8', quiet=TRUE)
# The supplied CRS is stored in .prj.txt rather than .prj. Read its actual WKT.
if (is.na(st_crs(v))) st_crs(v) <- paste(readLines(file.path(cache, 'Vegetacao_5000.prj.txt'), warn=FALSE), collapse='')
stopifnot(grepl('SIRGAS', st_crs(v)$wkt))
# Reconstruct the mapped original formation: surviving Mm/Ml polygons plus
# transformed polygons whose previous vegetation is pure Mixed Ombrophilous Forest.
# Mixed/contact classes (OM, EM, SM, NM) are intentionally not included.
previous <- !is.na(v$DSC_VEG_PR) & grepl('Mista', v$DSC_VEG_PR) & !grepl('/', v$DSC_VEG_PR)
keep <- previous | v$COD_CLASS_ %in% c('Mm','Ml')
forest <- st_transform(st_make_valid(v[keep,]), 4326)
# Use GADM 4.1, matching the administrative boundaries in the article.
# GADM data are obtained by each user and are not redistributed with the repository.
# Terms: https://gadm.org/license.html
boundaries <- file.path(map_cache,'figure_1_boundaries.rds')
if (!file.exists(boundaries)) {
  mun <- st_as_sf(geodata::gadm('BRA',level=2,path=map_cache,version='4.1'))
  states <- st_as_sf(geodata::gadm('BRA',level=1,path=map_cache,version='4.1'))
  mun <- mun[mun$NAME_1=='Santa Catarina',]
  names(mun)[names(mun)=='NAME_2'] <- 'name_muni'
  states$abbrev_state <- sub('^BR\\.','',states$HASC_1)
  mun <- st_transform(mun,4326); states <- st_transform(states,4326)
  sc <- states[states$NAME_1=='Santa Catarina',]
  sa <- rnaturalearth::ne_countries(continent='South America',returnclass='sf',scale='medium')
  saveRDS(list(mun=mun,states=states,sc=sc,sa=sa,src='gadm'),boundaries)
}
L <- readRDS(boundaries)
s <- read.csv(file.path(root,'data/sites.csv'),fileEncoding='UTF-8')
s$type <- factor(ifelse(s$foco,'Focal population','Comparison remnants'),levels=c('Focal population','Comparison remnants'))
pts <- st_as_sf(s,coords=c('lon','lat'),crs=4326,remove=FALSE)
joined <- st_join(st_transform(pts,31982),st_transform(L$mun[,'name_muni'],31982),join=st_within)
stopifnot(nrow(joined)==8, !anyNA(joined$name_muni))
names_mun <- unique(joined$name_muni)
selected_mun <- L$mun[L$mun$name_muni %in% names_mun,]
# Labels in geographic coordinates, with leaders to the unchanged site points.
labels <- data.frame(code=c('PB1','PB2','PAI','BJS','Vargem','Epagri','Alu','PARNAMUL'),
 x=c(-50.31,-49.96,-49.88,-49.99,-50.97,-50.48,-50.37,-50.54),
 y=c(-27.94,-27.79,-27.94,-28.40,-27.515,-27.705,-27.865,-27.805))
labels <- merge(s,labels,by='code',sort=FALSE)
labels$nome[labels$code=='PARNAMUL'] <- 'Parque Natural\nMunicipal'
mlab <- data.frame(x=c(-50.67083,-50.33,-50.08157,-49.66489),y=c(-27.60418,-28.11,-27.97967,-28.37415),
 label=c('S\u00c3O JOS\u00c9\nDO CERRITO','LAGES','PAINEL','BOM JARDIM\nDA SERRA'))
# Confirm each municipality label lies inside its named municipality.
mlab$expected <- c('S\u00e3o Jos\u00e9 do Cerrito','Lages','Painel','Bom Jardim da Serra')
mlab_sf <- st_as_sf(mlab,coords=c('x','y'),crs=4326,remove=FALSE)
mlab_check <- st_join(st_transform(mlab_sf,31982),st_transform(L$mun[,'name_muni'],31982))
stopifnot(all(mlab_check$expected==mlab_check$name_muni))
red <- '#BC4137'; blue <- '#286D9B'; green <- '#71A585'
main <- ggplot() +
 geom_sf(data=L$states,fill='#E5E5E5',color='#B5B5B5',linewidth=.25)+
 geom_sf(data=L$mun,fill='#FAFAF8',color='#D0D0CB',linewidth=.25)+
 geom_sf(data=selected_mun,fill='#F1F2ED',color='#9B9D96',linewidth=.4)+
 geom_sf(data=L$sc,fill=NA,color='#777B77',linewidth=.45)+
 geom_text(data=mlab,aes(x,y,label=label),size=3.0,fontface='italic',color='#727A73',lineheight=.95)+
 geom_segment(data=labels,aes(x=lon,y=lat,xend=x,yend=y),color='#8C918D',linewidth=.28)+
 geom_point(data=s,aes(lon,lat,shape=type,fill=type,size=type),color='white',stroke=.6)+
 geom_label(data=labels,aes(x,y,label=nome),size=3.05,fontface=ifelse(labels$foco,'bold','plain'),
  fill=scales::alpha('white',.88),color='#242A27',linewidth=0,label.padding=unit(.11,'lines'),lineheight=.96)+
 scale_shape_manual(values=c(24,21))+scale_fill_manual(values=c(red,blue))+
 scale_size_manual(values=c(4.4,3.2))+guides(size='none',shape=guide_legend(order=1,override.aes=list(size=c(4.4,3.2))),fill=guide_legend(order=1))+
 coord_sf(xlim=c(-51.13,-49.43),ylim=c(-28.49,-27.41),expand=FALSE)+
 annotation_scale(location='br',width_hint=.26,height=unit(.15,'cm'),text_cex=.65)+
 annotation_north_arrow(location='tl',which_north='true',height=unit(.85,'cm'),width=unit(.7,'cm'),style=north_arrow_fancy_orienteering())+
 labs(x=NULL,y=NULL)+theme_bw(base_size=10)+
 theme(panel.grid=element_line(color='#E5E7E2',linewidth=.2),legend.title=element_blank(),
  legend.position='bottom',legend.text=element_text(size=9),legend.margin=margin(0,0,0,0),
  plot.margin=margin(6,5,2,3),axis.text=element_text(size=8),panel.border=element_rect(color='#656A65',linewidth=.5))
br <- L$sa[L$sa$admin=='Brazil',]
locator <- ggplot()+geom_sf(data=L$sa,fill='#F0F0ED',color='#AEB2AB',linewidth=.2)+
 geom_sf(data=br,fill='#DDE1D8',color='#8A9186',linewidth=.3)+
 geom_sf(data=L$sc,fill=red,color=red,linewidth=.25)+
 annotate('text',x=-52,y=-10,label='BRAZIL',fontface='bold',size=3.1,color='#434B42')+
 coord_sf(xlim=c(-82,-34),ylim=c(-56,12),expand=FALSE)+theme_void()+
 theme(panel.border=element_rect(fill=NA,color='#8E958C',linewidth=.4),plot.margin=margin(2,2,2,2))
regional <- ggplot()+geom_sf(data=L$sa,fill='#F4F4F0',color='#B1B7AE',linewidth=.2)+
 geom_sf(data=L$states,fill='#F4F4F0',color='#ABB1A8',linewidth=.25)+
 geom_sf(data=forest,fill=green,color=NA)+geom_sf(data=L$sc,fill=NA,color='#555F51',linewidth=.5)+
 annotate('rect',xmin=-51.13,xmax=-49.43,ymin=-28.49,ymax=-27.41,fill=NA,color=red,linewidth=.65)+
 annotate('text',x=-53,y=-23.7,label='PR',size=3,color='#555F51')+
 annotate('text',x=-52.7,y=-27,label='SC',size=3,color='#555F51')+
 annotate('text',x=-53.5,y=-30.2,label='RS',size=3,color='#555F51')+
 coord_sf(xlim=c(-58,-40),ylim=c(-34,-19),expand=FALSE)+theme_void()+
 theme(panel.border=element_rect(fill=NA,color='#8E958C',linewidth=.4),plot.margin=margin(2,2,2,2))
key <- ggdraw()+draw_label('Original Araucaria Forest',x=.17,y=.77,hjust=0,size=9)+
 draw_label('IBGE vegetation map (1:5,000,000)',x=.05,y=.43,hjust=0,size=7.6,color='#5B6258')+
 draw_label('Red outline: main map extent',x=.05,y=.12,hjust=0,size=7.6,color='#5B6258')+
 draw_grob(grid::rectGrob(gp=grid::gpar(fill=green,col=NA)),x=.05,y=.66,width=.08,height=.22)
right <- plot_grid(locator,regional,key,ncol=1,rel_heights=c(1.1,1,.28))
fig <- plot_grid(main,right,nrow=1,rel_widths=c(2.42,1),align='h',axis='tb')
ggsave(file.path(out,'Figure_1.tiff'),fig,width=10,height=6.4,units='in',dpi=600,device=grDevices::tiff,type='cairo',compression='lzw',bg='white')
cat('Figure 1 exported as TIFF (600 dpi).\n')
