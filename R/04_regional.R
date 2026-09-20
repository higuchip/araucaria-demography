# Regional adult structure and Figure 6.
suppressPackageStartupMessages({ library(data.table) })
base <- normalizePath(getwd(), winslash="/")
ms <- fread(file.path(base, "data", "multisitio_araucaria.csv"))
ms[, BA := pi * (dap_cm/2)^2 / 10000]
per <- ms[, .(n = .N, BA = sum(BA)), by = .(sitio, ano)]
res <- per[, { o <- order(ano); yr <- ano[o]; n <- n[o]; BA <- BA[o]; T <- max(yr) - min(yr)
  list(Tyr = T, n1 = n[length(n)],
       dens = 100*((n[length(n)]/n[1])^(1/T) - 1),
       ba   = 100*((BA[length(BA)]/BA[1])^(1/T) - 1)) }, by = sitio]
nm <- c(PB1="Pedras Brancas 1", PB2="Pedras Brancas 2", PAI="Painel", BJS="Bom Jardim da Serra",
        Vargem="Vargem", Epagri="Amola Faca", Alu="Aluvial", PARNAMUL="PARNAMUL")
res[, nome := nm[sitio]]; res[, foco := sitio == "PB1"]
pal <- colorRampPalette(c("#bcd2e8","#08306b"))(100)
res[, ci := round(1 + 99*(Tyr - min(Tyr))/(max(Tyr) - min(Tyr)))]; res[, col := pal[ci]]
res[, cx := 1.0 + 2.0*sqrt(n1/max(n1))]

L <- list(
  pt = list(x="Mudan\u00e7a anual na densidade (%/ano)", y="Mudan\u00e7a anual na \u00e1rea basal (%/ano)",
            q1="AB \u2191 \u00b7 densidade \u2193\nmatura\u00e7\u00e3o / autodesbaste", q2="AB \u2191 \u00b7 densidade \u2191\nac\u00famulo / crescimento",
            q3="AB \u2193 \u00b7 densidade \u2193\nperda / substitui\u00e7\u00e3o", q4="AB \u2193 \u00b7 densidade \u2191\ningresso sem ac\u00famulo",
            mid="nenhum remanescente perdeu area basal",
            t1="Din\u00e2mica adulta de Araucaria angustifolia entre remanescentes protegidos",
            t2="(c\u00edrculo vermelho = \u00e1rea focal, com ciclo de vida completo)",
            ld="Dura\u00e7\u00e3o (anos)", ln="n\u00ba arauc\u00e1rias"),
  en = list(x="Annual change in density (%/yr)", y="Annual change in basal area (%/yr)",
            q1="BA \u2191 \u00b7 density \u2193\nmaturation / self-thinning", q2="BA \u2191 \u00b7 density \u2191\naccumulation / growth",
            q3="BA \u2193 \u00b7 density \u2193\nloss / replacement", q4="BA \u2193 \u00b7 density \u2191\nrecruitment without accumulation",
            mid="no remnant lost basal area",
            t1="Adult dynamics of Araucaria angustifolia across protected remnants",
            t2="(red circle = focal site, with full life cycle)",
            ld="Duration (years)", ln="No. of araucaria"))

draw <- function(k){
  par(mar = c(7.8,4.9,3,1.4), mgp = c(2.8,0.7,0))
  plot(res$dens, res$ba, xlim=c(-2.5,2.5), ylim=c(-6,6), type="n", asp=NA, xlab=k$x, ylab=k$y)
  abline(h=0, v=0, lty=2, col="grey60"); qc<-"grey48"; q<-.66
  text(-2.45, 5.8, k$q1, adj=c(0,1), cex=q, col=qc); text( 2.45, 5.8, k$q2, adj=c(1,1), cex=q, col=qc)
  text(-2.45,-5.8, k$q3, adj=c(0,0), cex=q, col=qc); text( 2.45,-5.8, k$q4, adj=c(1,0), cex=q, col=qc)
  text(0,-3.0, k$mid, col="grey55", cex=.82, font=3)
  points(res$dens, res$ba, pch=21, cex=res$cx, bg=res$col, col="white", lwd=1.4)
  with(res[foco==TRUE], points(dens, ba, pch=21, cex=cx, bg=NA, col="#C44E52", lwd=2.8))
  pos <- ifelse(res$foco,1,3); pos[res$sitio=="Vargem"]<-4; pos[res$sitio=="PARNAMUL"]<-2
  text(res$dens, res$ba, res$nome, pos=pos, cex=.78, font=ifelse(res$foco,2,1),
       col=ifelse(res$foco,"#C44E52","grey20"))
  mtext(k$t1, 3, 1.2, font=2, cex=.95); mtext(k$t2, 3, .2, cex=.78, col="grey40")
  Ts<-c(2,7,13); lc<-pal[round(1+99*(Ts-min(res$Tyr))/(max(res$Tyr)-min(res$Tyr)))]; ns<-c(25,90,155)
  par(xpd=NA); yL<-grconvertY(0.075,"ndc","user")
  legend(grconvertX(0.10,"ndc","user"), yL, horiz=TRUE, title=k$ld, legend=Ts, pt.bg=lc, pch=21, col="white", pt.cex=1.5, bty="n", cex=.78)
  legend(grconvertX(0.60,"ndc","user"), yL, horiz=TRUE, title=k$ln, legend=ns, pch=21, pt.bg="grey75", col="white", pt.cex=1.0+1.8*sqrt(ns/max(res$n1)), bty="n", cex=.78)
  par(xpd=FALSE)
}

dest <- file.path(base, 'out', 'figures')
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
tiff(file.path(dest,'Figure_6.tiff'),width=4364,height=4418,units='px',res=600,type='cairo',compression='lzw')
draw(L[['en']]); dev.off()
tables$regional_rates <- res[,.(sitio,nome,Tyr,n_final=n1,density_percent_year=dens,basal_area_percent_year=ba)]

sites <- fread(file.path(base,"data","sites.csv"))
regional <- merge(per,sites[,.(sitio=code,area_ha=parcelas*area_parcela_m2/10000)],by="sitio")
regional[,`:=`(density_ha=n/area_ha,basal_area_m2ha=BA/area_ha)]
tables$regional_structure <- regional
agg <- ms[,.(n_araucarias=uniqueN(individuo),first_census=min(ano),last_census=max(ano),inventories=uniqueN(ano)),by=sitio]
metadata <- merge(sites,agg,by.x="code",by.y="sitio")
metadata[,sampled_area_ha:=parcelas*area_parcela_m2/10000]
tables$Table_S1_1 <- metadata[match(c("PB1","PB2","Alu","PARNAMUL","Epagri","Vargem","PAI","BJS"),code)]
