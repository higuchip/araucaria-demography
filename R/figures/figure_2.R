src <- file.path('out')
dest <- file.path('out','figures')
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
mat <- as.data.frame(tables$matrix_elements)
death <- as.data.frame(tables$annual_death_probabilities)
mat <- subset(mat, model == 'integrated_uniform')
death <- subset(death, model == 'integrated_uniform')
stages <- c('C1', 'C2', 'C3', 'A')
stopifnot(nrow(mat) == 16L, nrow(death) == 4L)
rate <- function(from, to) {
  z <- mat[mat$from == from & mat$to == to, ]
  stopifnot(nrow(z) == 1L)
  z$median
}
pct <- function(x) sprintf(if (100*x < 1) '%.2f%%' else '%.1f%%', 100*x)

shown <- subset(mat, median > 0)
shown$kind <- ifelse(shown$from == 'A' & shown$to == 'C1', 'effective recruitment',
                     ifelse(shown$from == shown$to, 'stasis', 'transition'))
death$to <- 'Death'; death$kind <- 'mortality'
tab <- rbind(shown[,c('from','to','kind','median','lower','upper')],
             death[,c('from','to','kind','median','lower','upper')])
tab$unit <- ifelse(tab$kind == 'effective recruitment', 'seedlings adult^-1 year^-1', 'percent per year')
factor <- ifelse(tab$kind == 'effective recruitment', 1, 100)
tab[,c('median','lower','upper')] <- tab[,c('median','lower','upper')] * factor
tab$display <- mapply(function(m,k) if(k == 'effective recruitment') sprintf('%.1f',m) else pct(m/100), tab$median,tab$kind)
tab$shown_in_figure <- !(tab$from == 'C2' & tab$to == 'C1' | tab$from == 'C3' & tab$to == 'C2')
tables$figure_2_rates <- as.data.table(tab[tab$shown_in_figure,])

draw <- function() {
  par(mar=rep(0,4), family='sans', xaxs='i', yaxs='i', lend='round', ljoin='round')
  plot.new(); plot.window(xlim=c(0,180), ylim=c(25,112), asp=1)
  ink <- '#303A3B'; blue <- '#315F78'; light <- '#F0F4F2'; grey <- '#52605D'
  label <- function(x,y,s,size=9,bold=FALSE,col=ink) text(x,y,s,cex=size/12,font=if(bold) 2 else 1,col=col)
  curve_arrow <- function(p, col=ink, lty=1) {
    t <- seq(0,1,length.out=151)
    xy <- outer((1-t)^3,p[1,]) + outer(3*(1-t)^2*t,p[2,]) +
      outer(3*(1-t)*t^2,p[3,]) + outer(t^3,p[4,])
    lines(xy,lwd=1.25,col=col,lty=lty)
    arrows(xy[146,1],xy[146,2],xy[151,1],xy[151,2],length=.065,angle=25,lwd=1.25,col=col)
  }
  # The outer dashed arrow is recruitment, not a transition probability.
  curve_arrow(rbind(c(169,59),c(180,112),c(2,112),c(11,59)),col=blue,lty=2)
  label(90,109,'Effective recruitment',9.5,TRUE,blue)
  label(90,104,sprintf('%.1f seedlings per adult per year',rate('A','C1')),9,FALSE,blue)
  # C2 -> A connects annual census endpoints and can pass through C3 between censuses.
  curve_arrow(rbind(c(72,59),c(72,86),c(155,86),c(155,59)))
  label(113.5,84,pct(rate('C2','A')),9.5,TRUE)
  x <- c(18,66,114,162)
  # Retrogression remains in the model but is omitted from this simplified diagram.
  for(i in seq_along(x)) rect(x[i]-15,41,x[i]+15,59,col=light,border=grey,lwd=1.1)
  for(i in 1:4) label(x[i],54,c('C1','C2','C3','Adult')[i],11,TRUE)
  label(x[1],47,'0.20\u2013<1 m',8.5)
  label(x[2],47,'1\u20133 m',8.5)
  label(x[3],48,'> 3 m',8)
  label(x[3],43.8,'DBH < 5 cm',8)
  label(x[4],47,'DBH \u2265 5 cm',8.5)
  for(i in 1:3) {
    arrows(x[i]+15.7,50,x[i+1]-15.7,50,length=.065,angle=25,lwd=1.25,col=ink)
    label(mean(x[c(i,i+1)]),55,pct(rate(stages[i],stages[i+1])),9.5,TRUE)
  }
  # Stasis and mortality use the same convention in all four stages.
  # Marginal posterior medians need not sum to 100%.
  for(i in 1:4) {
    label(x[i],35,paste('Stasis',pct(rate(stages[i],stages[i]))),8)
    label(x[i],29,paste('Mortality',pct(death$median[death$from==stages[i]])),8)
  }
}
tiff(file.path(dest,'Figure_2.tiff'),width=180,height=87,units='mm',res=600,type='cairo',compression='lzw',pointsize=12)
draw(); dev.off()
