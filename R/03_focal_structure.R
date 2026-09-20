# Focal stage densities, adult structure, and Supplementary Section S2.
suppressPackageStartupMessages({ library(data.table) })
base <- normalizePath(getwd(), winslash="/")
dir.create(file.path(base, "out"), showWarnings = FALSE)

# ------------------------------------------------------------
# 1) ARBOREO (todas spp; 50 parcelas 10x20 = 200 m2 ~ 1 ha)
# ------------------------------------------------------------
ad <- fread(file.path(base, "data", "focal_arboreo.csv"))
ad[, ano := as.integer(ano)]
ad[, tree_id := paste(transecto, parcela, arvore, sep = "_")]
cat("==== ARBOREO (comunidade) ====\n")
cat("registros:", nrow(ad), " especies:", uniqueN(ad$especie),
    " arvores:", uniqueN(ad$tree_id),
    " censos:", paste(sort(unique(ad$ano)), collapse = "/"), "\n")
ad_ara <- ad[especie == "Araucaria angustifolia"]
cat("Araucaria adulta - arvores unicas:", uniqueN(ad_ara$tree_id), "\n")

# ------------------------------------------------------------
# 2) REGENERACAO (todas spp; classes 1/2/3; aninhada 5/10/20 m2; 2011-2014)
# ------------------------------------------------------------
reg <- fread(file.path(base, "data", "focal_regeneracao.csv"))
reg[, reg_id := paste(transecto, parcela, subparcela, individuo, sep = "_")]
reg[, ano := as.integer(ano)]
# d$rg: individuos (todas spp) com T,P,PC para mapear parcela <-> subparcela (PC)
rg <- unique(reg[, .(T = transecto, P = parcela, PC = subparcela, N = individuo, especie)])
# Araucaria -> formato longo (1 linha por individuo x ano; estado em {1,2,3,M,AD,NA})
rga <- reg[grepl("^Araucaria angustifolia", especie)]
reg_long <- rga[, .(reg_id, T = transecto, P = parcela, PC = subparcela,
                    ano, estado = classe)]
setorder(reg_long, reg_id, ano)
cat("\n==== REGENERACAO ====\n")
cat("especies:", uniqueN(reg$especie), " | Araucaria - individuos:", uniqueN(reg_long$reg_id),
    " | anos:", paste(sort(unique(reg_long$ano)), collapse = "/"), "\n")

# ------------------------------------------------------------
# 3) Densidades padronizadas (ind/ha) por classe e ano (Araucaria)
# ------------------------------------------------------------
n_parc <- 50
area_classe <- c("1" = 5, "2" = 10, "3" = 20)          # m2 por parcela
area_tot     <- area_classe * n_parc                    # m2 amostrados (total)
dens <- reg_long[estado %in% c("1","2","3"), .(n = .N), by = .(ano, classe = estado)]
dens[, dens_ha := n * 10000 / area_tot[classe]]
dens_wide <- dcast(dens, ano ~ classe, value.var = "dens_ha")
cat("\n==== DENSIDADE Araucaria (ind/ha) por classe/ano ====\n"); print(dens_wide)

ad_dens <- ad_ara[status %in% c("viva","recrutamento"),
                  .(adultos_ha = .N), by = ano][order(ano)]   # 50 parc x 200 m2 = 1 ha
cat("\n==== Adultos vivos Araucaria (ind/ha ~ contagem em 1 ha) ====\n"); print(ad_dens)

# ------------------------------------------------------------
# 4) Salvar
# ------------------------------------------------------------
tables$ad_ara <- copy(ad_ara)
tables$Table_1 <- merge(dens_wide,ad_dens,by="ano",all=TRUE)[order(ano)]
setnames(tables$Table_1,c("ano","1","2","3","adultos_ha"),c("year","C1_ind_ha","C2_ind_ha","C3_ind_ha","adults_ind_ha"))
species <- ad[ano==2008 & status=="viva", .(N=.N), by=especie][order(especie)]
stopifnot(sum(species$N)==1875L,nrow(species)==86L)
tables$Table_S1_2 <- species

ad[, dap_cm := as.numeric(dap_cm)]; ad[, ab_cm2 := as.numeric(ab_cm2)]
ad[, plot := paste(transecto,parcela)]
ar <- ad[especie=="Araucaria angustifolia"]
ha <- 1   # 50 x 200 m2 = 1 ha

traj <- function(x, lab){
  t <- x[status %in% c("viva","recrutamento"),
         .(N=.N, DAP_medio=mean(dap_cm,na.rm=TRUE),
           BA_m2ha=sum(ab_cm2,na.rm=TRUE)/10000/ha), by=ano][order(ano)]
  assign(if(lab=="ARAUCARIA") "adult_structure_araucaria" else "adult_structure_community",copy(t),envir=tables)
  cat("\n==== Trajetoria",lab,"(densidade, DAP medio, AB) ====\n"); print(t)
}
traj(ar, "ARAUCARIA")
traj(ad, "POVOAMENTO (todas spp)")

# ------------------------------------------------------------
# (2) Mortalidade assimetrica por tamanho - Araucaria
#   compara DAP inicial (2008) de quem MORREU ate 2021 vs SOBREVIVEU
# ------------------------------------------------------------
aw <- dcast(ar, plot+arvore ~ ano, value.var="status")
dapw <- dcast(ar, plot+arvore ~ ano, value.var="dap_cm")
viv08 <- aw[["2008"]]=="viva"; viv08[is.na(viv08)] <- FALSE
fate <- data.table(dap08=dapw[["2008"]][viv08],
                   morreu = aw[["2021"]][viv08]=="morta")
fate <- fate[!is.na(dap08)]
cat("\n==== Mortalidade por tamanho (Araucaria, coorte 2008->2021) ====\n")
cat(sprintf("DAP medio: sobreviventes=%.1f cm (n=%d) | mortos=%.1f cm (n=%d)\n",
            mean(fate[morreu==FALSE]$dap08), nrow(fate[morreu==FALSE]),
            mean(fate[morreu==TRUE]$dap08),  nrow(fate[morreu==TRUE])))
print(suppressWarnings(wilcox.test(dap08 ~ morreu, data=fate)))

# ------------------------------------------------------------
# (3) Mortalidade denso-dependente entre parcelas (Araucaria)
#   taxa de mortos 2008-2021 ~ densidade/AB conspecifica inicial
# ------------------------------------------------------------
pl <- ar[ano==2008 & status=="viva",
         .(N0=.N, BA0=sum(ab_cm2,na.rm=TRUE)), by=plot]
mort <- ar[status=="morta", .(mortos=uniqueN(arvore)), by=plot]
pl <- merge(pl, mort, by="plot", all.x=TRUE); pl[is.na(mortos),mortos:=0]
pl[, taxa := mortos/N0]
cat("\n==== Mortalidade denso-dependente (Araucaria, ", nrow(pl)," parcelas) ====\n", sep="")
cat(sprintf("Spearman taxa_mort ~ densidade conspecifica: rho=%+.3f p=%.3f\n",
            cor(pl$N0,pl$taxa,method="spearman"),
            suppressWarnings(cor.test(pl$N0,pl$taxa,method="spearman"))$p.value))
cat(sprintf("Spearman taxa_mort ~ AB conspecifica:        rho=%+.3f p=%.3f\n",
            cor(pl$BA0,pl$taxa,method="spearman"),
            suppressWarnings(cor.test(pl$BA0,pl$taxa,method="spearman"))$p.value))

# ------------------------------------------------------------
# (4) Spatial density-size relationship: community plots in 2008
#   An among-plot relationship, not a temporal self-thinning trajectory.
# ------------------------------------------------------------
st <- ad[ano==2008 & status=="viva",
         .(N=.N, DAPm=mean(dap_cm,na.rm=TRUE)), by=plot][N>=3]
fit <- lm(log(N) ~ log(DAPm), data=st)
cat("\n==== Linha de autodesbaste (povoamento, ", nrow(st)," parcelas) ====\n", sep="")
cat(sprintf("log(N) ~ log(DAP medio): inclinacao = %.2f  p=%.3f  R2=%.2f\n",
            coef(fit)[2], summary(fit)$coefficients[2,4], summary(fit)$r.squared))
cat("Spatial association; not a temporal self-thinning slope.\n")

fate[,dbh_class:=cut(dap08,c(5,15,25,35,Inf),right=FALSE)]
tables$Table_S2_2 <- fate[,.(n=.N,deaths=sum(morreu),mortality_percent=100*mean(morreu)),by=dbh_class]
tables$plot_mortality <- pl
tables$Table_S2_4 <- data.table(slope=coef(fit)[2],p=summary(fit)$coefficients[2,4],R2=summary(fit)$r.squared)

ar <- ad[especie=="Araucaria angustifolia"]
ar[, dap_cm := as.numeric(dap_cm)]; ar[, id := paste(transecto,parcela,arvore)]

brk <- c(5,15,25,35,45,Inf); lab <- c("5-15","15-25","25-35","35-45","45+")
ar[, cls := cut(dap_cm, brk, lab, right=FALSE)]

# (1) distribuicao diametrica por ano (vivos) -------------------------
dd <- dcast(ar[status %in% c("viva","recrutamento")], cls ~ ano, fun.aggregate=length)
cat("==== Distribuicao diametrica (n vivos por classe/ano) ====\n"); print(dd)

# (2) crescimento e mortalidade anual por classe ----------------------
ar2 <- ar[order(id,ano)]
ar2[, dap_next := shift(dap_cm,-1), by=id]
ar2[, ano_next := shift(ano,-1), by=id]
ar2[, st_next := shift(status,-1), by=id]
ar2[, dt := ano_next-ano]
seg <- ar2[!is.na(ano_next) & status=="viva"]
cat("\n==== Crescimento DAP (cm/ano) e mortalidade anual por classe ====\n")
gm <- seg[, .(n=.N,
              cresc_cm_ano = mean((dap_next-dap_cm)/dt, na.rm=TRUE),
              mort_anual = 1 - (sum(st_next=="viva")/.N)^(1/mean(dt))), by=cls][order(cls)]
print(gm)


tables$adult_diameter_classes <- dd
tables$adult_diameter_growth_mortality <- gm

tables$Table_S2_1 <- merge(tables$adult_structure_araucaria,tables$adult_structure_community,
 by="ano",suffixes=c("_araucaria","_community"))
tables$Table_S2_3 <- rbindlist(lapply(c("N0","BA0"),function(v){
 test <- suppressWarnings(cor.test(pl[[v]],pl$taxa,method="spearman"))
 data.table(predictor=if(v=="N0") "Initial conspecific density" else "Initial conspecific basal area",
 rho=unname(test$estimate),p=test$p.value,n_plots=nrow(pl))
}))
w <- suppressWarnings(wilcox.test(dap08 ~ morreu,data=fate))
tables$focal_text <- data.table(metric=c("mean_initial_DBH_survivors","mean_initial_DBH_deaths","size_mortality_W","size_mortality_p"),
 value=c(mean(fate[morreu==FALSE]$dap08),mean(fate[morreu==TRUE]$dap08),unname(w$statistic),w$p.value),
 unit=c("cm","cm","statistic","probability"))
