# Run from the repository root: Rscript R/run_all.R
packages <- c("data.table", "popbio", "Rage", "ggplot2", "sf", "ggspatial", "cowplot", "geodata", "rnaturalearth", "rnaturalearthdata")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install these R packages first: ",paste(missing,collapse=", "))
stopifnot(file.exists("data/focal_arboreo.csv"))
# Intermediate objects remain in memory and are shared with the plotting scripts.
tables <- new.env(parent=emptyenv())
scripts <- c("R/01_annual_model.R","R/03_focal_structure.R","R/04_regional.R",
             sprintf("R/figures/figure_%d.R",1:5),"R/figures/figure_S1.R")
for(script in scripts){
 cat("\nRunning ",script,"\n",sep="")
 module <- new.env(parent=globalenv())
 module$tables <- tables
 source(script,local=module,encoding="UTF-8")
}

# Tables follow the numbering in the manuscript and Supplementary Material.
s <- tables$comparison_summary
metrics <- c("lambda","R0","T_R0","adult_expected_years","C1_residence","coast50","coast10","damping")
tab2 <- data.table::copy(s[match(metrics,metric),.(metric,median,lower,upper)])
tab2[,parameter:=c("Asymptotic annual population growth rate","Net reproductive rate",
 "Generation time","Expected adult residence","Expected cumulative residence in C1",
 "Time to 50% without adult ingrowth","Time to 10% without adult ingrowth","Damping ratio")]
tab2[,unit:=c("ratio","ratio",rep("years",5),"ratio")]
amp <- tables$transient_summary[metric=="maxamp_0_400",.(metric,median,lower,upper)]
amp[,`:=`(parameter="Maximum transient amplification (2012 structure, 0-400 years)",unit="ratio")]
tables$Table_2 <- data.table::rbindlist(list(tab2,amp))
data.table::setcolorder(tables$Table_2,c("parameter","unit","median","lower","upper","metric"))
output <- file.path("out","tables")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
for(name in c("Table_1","Table_2","Table_S1_1","Table_S1_2","Table_S2_1","Table_S2_2","Table_S2_3","Table_S2_4")){
 data.table::fwrite(tables[[name]],file.path(output,paste0(name,".csv")),na="NA")
}

cat("\nFinished: eight CSV tables and seven TIFF figures (600 dpi) in out/.\n")
invisible(tables)
