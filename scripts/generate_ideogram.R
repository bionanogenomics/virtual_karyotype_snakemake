args <- commandArgs(trailingOnly=TRUE)
print(args)
cyto <- args[1]
genome <- args[2]
orientation <- args[3]
kprect <- args[4]
processed <- args[5]
sample_handle <- args[6]

library(karyoploteR)


custom.cytobands <- toGRanges(cyto)
custom.genome <- toGRanges(genome)
orientation_frame = read.table(orientation, sep='\t', header=T)
kprect_frame = read.table(kprect, sep='\t', header=T)


mapColors <- function(
   contig_list, brewerPaletteName= 'Set3', reverse=FALSE, colors=NULL
) {
   if (! is.null(colors)) {
      labels <- colors
   } else {
      labels= colorRampPalette(RColorBrewer::brewer.pal( length(unique(contig_list)), brewerPaletteName))(length(unique(contig_list)))
      names(labels) <- unique(contig_list)
   }
   if (reverse) {
      labels <- rev(labels)
   }
   contig_list<-as.character(contig_list)
   labelCol <- unlist(lapply(contig_list, FUN = function(x){labels[x]}))
   return(labelCol);
}

plot_kprects <- function(
   sub_kprect, path, color_mapper, kp
   ){
      print(sub_kprect)
      x0 <- as.numeric(sub_kprect[2])
      x1 <- as.numeric(sub_kprect[3])
      chrom <- sub_kprect[9]
      midpoint <- as.numeric(sub_kprect[20])
      kpRect(kp, chr=path, x0=x0, x1=x1, y0=0.0, y1=1.0, col=color_mapper[chrom], lwd=0.1)
      kpText(kp, chr=path, x=midpoint, y=0.5, labels=chrom, cex=1.25, srt=90)
}

plot_orientation <- function(
   sub_orientation, path, kp
){
   x0 = as.numeric(sub_orientation[2])
   x1 = as.numeric(sub_orientation[3])
   col = sub_orientation[20]
   kpArrows(kp, chr=path, x0=x0, x1=x1, y0=0.5, y1=0.5, data.panel=2, lwd=1, length=0.05, angle=20, col=col)
}

plot_ideogram <- function(
   path, fig_out, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper 
){
   pdf(fig_out)
   kp <- plotKaryotype(chromosomes=path, genome=sub_genome, plot.type=2, cytobands=sub_cytoband, plot.params=pp)
   kpAddBaseNumbers(kp)
   kpAddCytobandLabels(kp, force.all=TRUE, srt=90, col='#2C02FD', cex=0.5)
   apply(sub_kprect, 1, plot_kprects, path=path, color_mapper=color_mapper, kp=kp)
   apply(sub_orientation, 1, plot_orientation, path=path, kp=kp)
   dev.off()
}

plot_total_ideogram <- function(
   path, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp
){
   apply(sub_kprect, 1, plot_kprects, path=path, color_mapper=color_mapper, kp=kp)
   apply(sub_orientation, 1, plot_orientation, path=path, kp=kp)
}

chunk.2 <- function(x, n, force.number.of.groups = TRUE, len = length(x), groups = trunc(len/n), overflow = len%%n) { 
  if(force.number.of.groups) {
    f1 <- as.character(sort(rep(1:n, groups)))
    f <- as.character(c(f1, rep(n, overflow)))
  } else {
    f1 <- as.character(sort(rep(1:groups, n)))
    f <- as.character(c(f1, rep("overflow", overflow)))
  }
  
  g <- split(x, f)
  
  if(force.number.of.groups) {
    g.names <- names(g)
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
  } else {
    g.names <- names(g[-length(g)])
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
    g.names.ordered <- c(g.names.ordered, "overflow")
  }
  
  return(g[g.names.ordered])
}


pp <- getDefaultPlotParams(plot.type=2)
pp$leftmargin <- 0.3
pp$data2height <- 30

color_mapper <- mapColors(unique(custom.cytobands$chrom))
genome_split = split(custom.genome,custom.genome@seqnames)
cytoband_split = split(custom.cytobands, custom.cytobands@seqnames)
orientation_split = split(orientation_frame, orientation_frame$chr)
kprect_split = split(kprect_frame, kprect_frame$chr)
for(i in names(genome_split)){
   print(i)
   path_handle = strsplit(i, " ")[[1]][1]
   fig_out = paste(sample_handle, path_handle, '.pdf',sep='')
   sub_cytoband <- custom.cytobands[custom.cytobands@seqnames==i]
   sub_genome <- custom.genome[custom.genome@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_ideogram(i, fig_out, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper)
}

################
pp <- getDefaultPlotParams(plot.type=2)
pp$leftmargin <- 0.3
pp$data2height <- 100
pp$ideogramheight <- 550 
pp$dataideogrammin <- -.75
pp$dataideogrammax <- .75
pp$data2inmargin <- 20
pp$data1inmargin <- 20
pp$data1height <- 300
pp$tick.len <- 30
pp$data1outmargin <- 20
pp$data2outmargin <- 20
pp$data1max=1
pp$data1min=0
pp$bottommargin = 30
pp$topmargin = 30

cyto <- data.frame(custom.cytobands)
sub_cyto <- cyto[c('seqnames','chrom')]
non_dup <- sub_cyto[!duplicated(sub_cyto),]
chromosomes <- non_dup$seqnames[!duplicated(non_dup$seqnames)]
chrom_order <- rev(as.character(chromosomes))
chrom_split_order <- chunk.2(chrom_order, 2, force.number.of.groups=T)
seqlevels(custom.genome) <- chrom_order
custom.genome <- sort(custom.genome)
seqlevels(custom.cytobands,) <- chrom_order
levels(custom.genome@seqnames) <- chrom_order
levels(custom.cytobands@seqnames) <- chrom_order

contig_split <- round(length(chrom_order)/2)

chrom_split_order <- chunk.2(chrom_order, 2, force.number.of.groups=T)

genome_part_1 <- custom.genome[seqnames(custom.genome) %in% chrom_split_order$`2`]
cytoband_part_1 <- subsetByOverlaps(custom.cytobands,genome_part_1)
seqlevels(genome_part_1) <- seqlevelsInUse(genome_part_1)
levels(genome_part_1) <- seqlevels(genome_part_1)
seqlevels(cytoband_part_1) <- seqlevelsInUse(cytoband_part_1)
levels(cytoband_part_1) <- seqlevels(cytoband_part_1)

fig_out = paste(sample_handle, 'karyotype_split_1_of_2', '.pdf',sep='')
pdf(fig_out, width=10,height=20,pointsize=1)
genome_split = split(genome_part_1,genome_part_1@seqnames)
kp <- plotKaryotype(chromosomes=chrom_split_order$`2`, genome=genome_part_1, plot.type=2, cytobands=cytoband_part_1, plot.params=pp, lwd=0.1, cex=3.5)
kpAddBaseNumbers(kp,cex=.8)
kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
for(i in names(genome_split)){
   sub_cytoband <- cytoband_part_1[cytoband_part_1@seqnames==i]
   sub_genome <- genome_part_1[genome_part_1@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
}
dev.off()


genome_part_2 <- custom.genome[seqnames(custom.genome) %in% chrom_split_order$`1`]
cytoband_part_2 <- subsetByOverlaps(custom.cytobands,genome_part_2)
seqlevels(genome_part_2) <- seqlevelsInUse(genome_part_2)
levels(genome_part_2) <- seqlevels(genome_part_2)
seqlevels(cytoband_part_2) <- seqlevelsInUse(cytoband_part_2)
levels(cytoband_part_2) <- seqlevels(cytoband_part_2)

fig_out = paste(sample_handle, 'karyotype_split_2_of_2', '.pdf',sep='')
pdf(fig_out, width=10,height=20,pointsize=1)
genome_split = split(genome_part_2,genome_part_2@seqnames)
kp <- plotKaryotype(chromosomes=chrom_split_order$`1`, genome=genome_part_2, plot.type=2, cytobands=cytoband_part_2, plot.params=pp, cex=3.5)
kpAddBaseNumbers(kp,cex=.8)
kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
for(i in names(genome_split)){
   sub_cytoband <- cytoband_part_2[cytoband_part_2@seqnames==i]
   sub_genome <- genome_part_2[genome_part_2@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
}
dev.off()


x <- data.frame()
write.table(x, file=processed,col.names=FALSE)
