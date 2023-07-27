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


####### TEST
loopAndApply <- function(valueList, genome, plot.params) {
  results <- vector("numeric", length(valueList))  # Initialize an empty vector for results
  
  for (i in seq_along(valueList)) {
    value <- valueList[i]
    result <- ideoMid(value, genome, plot.params)  # Call ideoMid function with value, genome, and plot.params
    results[i] <- result  # Store the result in the vector
  }
  
  # Create a data frame with the results
  result_df <- data.frame(value = valueList, result = results)
  
  return(result_df)
}

#loopAndApply(kp$chromosomes,genome, kp$plot.params)
ideoMid <- function(chr, genome, plot.params) {
  pp <- plot.params
  chr.height <- getChrHeight_2HorizDataAboveAndBelowIdeogram(pp)
  chr.names <- GenomeInfoDb::seqlevels(genome)
  chrs <- c(length(chr.names):1)
  names(chrs) <- chr.names
  chr.num <- chrs[chr]
  return(pp$bottommargin + (chr.num - 1) * chr.height +
         pp$data2outmargin + pp$data2height + pp$data2inmargin + pp$ideogramheight/2)
}

getChrHeight_2HorizDataAboveAndBelowIdeogram <- function(pp) {
  chr.height <- (pp$data2outmargin + pp$data2height + pp$data2inmargin
                 + pp$ideogramheight
                 + pp$data1inmargin + pp$data1height + pp$data1outmargin)
  return(chr.height)
}


# Function to summarize data based on unique values in the 'chrom' column
# summarize_data_by_chrom <- function(coord_frame, pp) {
#   # Find unique values in the 'chrom' column
#   unique_chrom_values <- unique(coord_frame$chrom)
  
#   # Initialize an empty data frame to store the summarized results
#   summary_df <- data.frame(chrom = character(),
#                            max_x = numeric(),
#                            avg_y_mid_norm_adj = numeric(),
#                            stringsAsFactors = FALSE)
#   max_x <- max(coord_frame$x) + .04
#   # Iterate through each unique 'chrom' value
#   for (chrom_val in unique_chrom_values) {
#     # Subset the data frame for the current 'chrom' value
#     subset_df <- coord_frame[coord_frame$chrom == chrom_val, ]
    
#     # Set the 'max_x' value for all rows in the current 'chrom' group to the calculated maximum
#     subset_df$max_x <- max_x
    
#     # Calculate the average value of 'y_mid_norm_adj'
#     avg_y_mid_norm_adj <- mean(subset_df$y_mid_norm_adj)
    
#     # Create a new data frame with the 'chrom', 'max_x', and 'avg_y_mid_norm_adj' values for the current 'chrom' value
#     current_summary <- data.frame(chrom = chrom_val,
#                                   max_x = max_x,
#                                   avg_y_mid_norm_adj = avg_y_mid_norm_adj)
    
#     # Append the current_summary to the summary_df
#     summary_df <- rbind(summary_df, current_summary)
#   }
#   summary_df$avg_y_mid_norm_adj <- calculate_values_df_updated(coord_frame, pp)
  
#   return(summary_df)
# }


plot_kprects <- function(
   sub_kprect, path, color_mapper, kp
   ){
      x0 <- as.numeric(sub_kprect[2])
      x1 <- as.numeric(sub_kprect[3])
      chrom <- sub_kprect[9]
      midpoint <- as.numeric(sub_kprect[20])
      kpRect(kp, chr=path, x0=x0, x1=x1, y0=0.0, y1=.5, col=color_mapper[chrom], lwd=0.1)
      kpText(kp, chr=path, x=midpoint, y=0.25, labels=chrom, cex=2, srt=90)
}

merge_frames <- function(df1, df2) {
  # Perform the merge based on the row names of df1
  merged_df <- merge(df1, df2, by = "row.names", all.x = TRUE, sort=FALSE)
  
  # Set the row names to the first column and remove the extra 'Row.names' column
  rownames(merged_df) <- merged_df$Row.names
  merged_df$Row.names <- NULL
  
  return(merged_df)
}

get_coords <- function(
   kp, cyto_first
   ){
      chr.labels <- kp$chromosomes
      ccf <- kp$coord.change.function
      coords = ccf(chr=chr.labels,x=end(kp$plot.region), y=start(kp$plot.region),data.panel=2)
      coord_frame = data.frame(coords)
      coord_frame$y_norm = (coord_frame$y / max(coord_frame$y))
      coord_frame$x_norm = coord_frame$x + 0.025
      genome = kp$plot.region
      mid_ideogram <- loopAndApply(kp$chromosomes, genome, kp$plot.params)
      coord_frame$y_mid <-mid_ideogram$result
      coord_frame$y_mid_norm = (coord_frame$y_mid / max(coord_frame$y_mid))
      coord_frame$y_mid_norm_adj = coord_frame$y_mid_norm + (coord_frame$y_mid_norm - coord_frame$y_norm)
      coord_frame <- merge_frames(coord_frame, cyto_first)
      return(coord_frame)
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

####### split into parts
get_first_rows_per_seqnames <- function(dataframe) {
  # Find unique values in the 'seqnames' column
  unique_seqnames <- unique(dataframe$seqnames)
  
  # Initialize an empty list to store the subset dataframes
  subset_dfs <- list()
  
  # Iterate through each unique 'seqnames' value
  for (seqname in unique_seqnames) {
    # Subset the dataframe for the current 'seqnames' value and take the first row
    subset_df <- dataframe[dataframe$seqnames == seqname, ][1, ]
    
    # Append the subset dataframe to the list
    subset_dfs[[seqname]] <- subset_df
  }
  
  # Combine the subset dataframes into a single dataframe
  result_df <- do.call(rbind, subset_dfs)
  
  return(result_df)
}

plot_labels_by_chrom <- function(data_frame) {
  # Get the number of rows in the data frame
  num_rows <- nrow(data_frame)
  
  # Create individual scatter plots for each row
  for (i in 1:num_rows) {
    # Extract data for the current row
    chrom_value <- data_frame$chrom[i]
    max_x_value <- data_frame$max_x[i]
    yadj <- data_frame$yadj[i]
    
    # Create a scatter plot with text labels for the current row
    text(x=max_x_value, y=yadj, labels = chrom_value, cex=5, srt=90, pos = 4, col = "blue")
  }
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



cyto_first <- get_first_rows_per_seqnames(sub_cyto)

part1_5 <- droplevels(unique(cyto_first[cyto_first$chrom %in% c(1,2,3,4,5),]$seqnames))
part6_12 <- droplevels(unique(cyto_first[cyto_first$chrom %in% c(6,7,8,9,10,11,12),]$seqnames))
part13_18 <- droplevels(unique(cyto_first[cyto_first$chrom %in% c(13,14,15,16,17,18),]$seqnames))
part19_Y <- droplevels(unique(cyto_first[cyto_first$chrom %in% c(19,20,21,22,'X','Y'),]$seqnames))


genome_part1_5 <- custom.genome[seqnames(custom.genome) %in% part1_5]
cytoband_part1_5 <- subsetByOverlaps(custom.cytobands,genome_part1_5)
seqlevels(genome_part1_5) <- seqlevelsInUse(genome_part1_5)
levels(genome_part1_5) <- seqlevels(genome_part1_5)
seqlevels(cytoband_part1_5) <- seqlevelsInUse(cytoband_part1_5)
levels(cytoband_part1_5) <- seqlevels(cytoband_part1_5)

genome_part6_12 <- custom.genome[seqnames(custom.genome) %in% part6_12]
cytoband_part6_12 <- subsetByOverlaps(custom.cytobands,genome_part6_12)
seqlevels(genome_part6_12) <- seqlevelsInUse(genome_part6_12)
levels(genome_part6_12) <- seqlevels(genome_part6_12)
seqlevels(cytoband_part6_12) <- seqlevelsInUse(cytoband_part6_12)
levels(cytoband_part6_12) <- seqlevels(cytoband_part6_12)

genome_part13_18 <- custom.genome[seqnames(custom.genome) %in% part13_18]
cytoband_part13_18 <- subsetByOverlaps(custom.cytobands,genome_part13_18)
seqlevels(genome_part13_18) <- seqlevelsInUse(genome_part13_18)
levels(genome_part13_18) <- seqlevels(genome_part13_18)
seqlevels(cytoband_part13_18) <- seqlevelsInUse(cytoband_part13_18)
levels(cytoband_part13_18) <- seqlevels(cytoband_part13_18)

genome_part19_Y <- custom.genome[seqnames(custom.genome) %in% part19_Y]
cytoband_part19_Y <- subsetByOverlaps(custom.cytobands,genome_part19_Y)
seqlevels(genome_part19_Y) <- seqlevelsInUse(genome_part19_Y)
levels(genome_part19_Y) <- seqlevels(genome_part19_Y)
seqlevels(cytoband_part19_Y) <- seqlevelsInUse(cytoband_part19_Y)
levels(cytoband_part19_Y) <- seqlevels(cytoband_part19_Y)

# Function to calculate descending values from 1.0 (exclusive) based on the number of unique entries in 'chrom'
calculate_descending_values <- function(coord_frame) {
  # Find unique values in the 'chrom' column
  unique_chrom_values <- unique(coord_frame$chrom)

  # Calculate the number of unique 'chrom' values
  num_unique_chrom_values <- length(unique_chrom_values)

  # Calculate the increment to ensure values won't be 0 or 1
  increment <- 1 / (num_unique_chrom_values + 1)

  # Generate descending values from 1.0 (exclusive) based on the number of unique 'chrom' values
  descending_values <- seq(1 - increment, increment, length.out = num_unique_chrom_values)

  return(descending_values)
}

# Function to calculate values and return a data frame
calculate_values_df <- function(coord_frame) {
  # Find unique values in the 'chrom' column
  unique_chrom_values <- unique(coord_frame$chrom)

  # Calculate the total number of entries in 'chrom'
  total_entries <- nrow(coord_frame)

  # Create an empty data frame to store the results
  result_df <- data.frame(chrom = character(),
                          value = numeric(),
                          stringsAsFactors = FALSE)

  # Iterate through each unique value in 'chrom'
  n_mult = length(unique_chrom_values)
  for (chrom_val in unique_chrom_values) {
    n_mult <- n_mult - 1
    # Calculate the number of occurrences of the current chrom value
    num_occurrences <- sum(coord_frame$chrom == chrom_val)

    # Calculate the value based on the number of occurrences and total number of entries
    value <- num_occurrences / total_entries
    
    # Create a new row for the result data frame
    result_row <- data.frame(chrom = chrom_val,
                             value = ((value * n_mult) + (value/num_occurrences)),
                             stringsAsFactors = FALSE)

    # Add the result row to the result data frame
    result_df <- rbind(result_df, result_row)
  }

  return(result_df$value)
}

group_continuous <- function(df, group_column) {
  # Generate helper vector that increments when group_column value changes
  helper <- c(0, cumsum(as.numeric(diff(df[[group_column]])) != 0))

  # Split the data frame based on the helper vector
  split_df <- split(df, helper)

  return(split_df)
}

# Function to calculate values and return a data frame
calculate_coordinate_space <- function(coord_frame, pp) {
      # Find unique values in the 'chrom' column
   unique_chrom_values <- unlist(as.list(as.character(coord_frame$chrom)))
   topmargin <- pp$topmargin
   bottommargin <- pp$bottommargin
   chrom_width <- calc_chrom_width(pp)

   # Calculate the total number of entries in 'chrom'
   total_entries <- nrow(coord_frame)

      # Create an empty data frame to store the results
   result_df <- data.frame(chrom = character(),
                           yadj = numeric(),
                           max_x = numeric(),
                           stringsAsFactors = FALSE)

   # Iterate through each unique value in 'chrom'
   n_mult = length(unique_chrom_values)
   n_total = length(unique_chrom_values)
   total_width <- (chrom_width * n_total)
   total_plot <- topmargin + bottommargin + total_width
   off_set <- 0 
   split_frames <- group_continuous(df=coord_frame,group_column='chrom')
   max_x <- max(coord_frame$x) + .035
   for (frame in split_frames) {
      if(n_mult == n_total){
      top_offset = topmargin # accounts for the top margin in first calculation
      } else{
      top_offset = 0
      }
      # Calculate the number of occurrences of the current chrom value
      num_occurrences <- dim(frame)[[1]]

      ind_chrom_width <- num_occurrences * chrom_width
      middle_chrom_width <- ind_chrom_width/2
      

      # Create a new row for the result data frame
      result_row <- data.frame(chrom = unique(frame$chrom),
                              yadj = 1 - (middle_chrom_width + off_set)/total_plot,
                              max_x = max_x,
                              stringsAsFactors = FALSE)
      off_set = off_set + ind_chrom_width + top_offset

      # Add the result row to the result data frame
      result_df <- rbind(result_df, result_row)
   }

   return(result_df)
}


# Function to calculate values and return a data frame
calc_chrom_width <- function(pp) {
  # Below needs to be calculated for number of paths in coord frame 
  data1outmargin <- pp$data1outmargin 
  data2outmargin <- pp$data2outmargin 
  data1inmargin <- pp$data1inmargin
  data2inmargin <- pp$data2inmargin 
  ideogramheight <- pp$ideogramheight
  data1height <- pp$data1height
  data2height <- pp$data1height

  chrom_width <- data1outmargin + data2outmargin + data1inmargin + data2inmargin + ideogramheight + data1height + data2height
  return(chrom_width)
}

add_labels <- function(kp,cyto_first,pp) {
   coord_frame <- get_coords(kp,cyto_first)
   coord_frame$chrom <- droplevels(coord_frame$chrom)
   chrom_coords<-calculate_coordinate_space(coord_frame,pp)
   plot_labels_by_chrom(chrom_coords)   
}

fig_out = paste(sample_handle, 'karyotype_split_1_of_4', '.pdf',sep='')
pdf(fig_out, width=15,height=40,pointsize=1)
genome_split = split(genome_part1_5,genome_part1_5@seqnames)
kp <- plotKaryotype(chromosomes=levels(genome_part1_5@seqnames), genome=genome_part1_5, plot.type=2, cytobands=cytoband_part1_5, plot.params=pp, lwd=0.1, cex=3.5)
kpAddBaseNumbers(kp,cex=.8)
kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
for(i in names(genome_split)){
   sub_cytoband <- cytoband_part1_5[cytoband_part1_5@seqnames==i]
   sub_genome <- genome_part1_5[genome_part1_5@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
}
add_labels(kp, cyto_first, pp)
dev.off()

fig_out = paste(sample_handle, 'karyotype_split_2_of_4', '.pdf',sep='')
pdf(fig_out, width=15,height=40,pointsize=1)
genome_split = split(genome_part6_12,genome_part6_12@seqnames)
kp <- plotKaryotype(chromosomes=levels(genome_part6_12@seqnames), genome=genome_part6_12, plot.type=2, cytobands=cytoband_part6_12, plot.params=pp, lwd=0.1, cex=3.5)
kpAddBaseNumbers(kp,cex=.8)
kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
for(i in names(genome_split)){
   sub_cytoband <- cytoband_part6_12[cytoband_part6_12@seqnames==i]
   sub_genome <- genome_part6_12[genome_part6_12@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
}
add_labels(kp, cyto_first, pp)
dev.off()

fig_out = paste(sample_handle, 'karyotype_split_3_of_4', '.pdf',sep='')
pdf(fig_out, width=15,height=40,pointsize=1)
genome_split = split(genome_part13_18,genome_part13_18@seqnames)
kp <- plotKaryotype(chromosomes=levels(genome_part13_18@seqnames), genome=genome_part13_18, plot.type=2, cytobands=cytoband_part13_18, plot.params=pp, lwd=0.1, cex=3.5)
kpAddBaseNumbers(kp,cex=.8)
kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
for(i in names(genome_split)){
   sub_cytoband <- cytoband_part13_18[cytoband_part13_18@seqnames==i]
   sub_genome <- genome_part13_18[genome_part13_18@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
}
add_labels(kp, cyto_first, pp)
dev.off()

fig_out = paste(sample_handle, 'karyotype_split_4_of_4', '.pdf',sep='')
pdf(fig_out, width=15,height=40,pointsize=1)
genome_split = split(genome_part19_Y,genome_part19_Y@seqnames)
kp <- plotKaryotype(chromosomes=levels(genome_part19_Y@seqnames), genome=genome_part19_Y, plot.type=2, cytobands=cytoband_part19_Y, plot.params=pp, lwd=0.1, cex=3.5)
kpAddBaseNumbers(kp,cex=.8)
kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
for(i in names(genome_split)){
   sub_cytoband <- cytoband_part19_Y[cytoband_part19_Y@seqnames==i]
   sub_genome <- genome_part19_Y[genome_part19_Y@seqnames==i] 
   sub_orientation = orientation_split[i][[1]]
   sub_kprect = kprect_split[i][[1]]
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
}
add_labels(kp, cyto_first, pp)
dev.off()


# ###################
# fig_out = paste('test_final_applied_v5_', 'karyotype_split_1_of_2', '.pdf',sep='')
# pdf(fig_out, width=15,height=40,pointsize=1)
# genome_split = split(genome_part_1,genome_part_1@seqnames)
# kp <- plotKaryotype(chromosomes=chrom_split_order$`2`, genome=genome_part_1, plot.type=2, cytobands=cytoband_part_1, plot.params=pp, lwd=0.1, cex=3.5)
# kpAddBaseNumbers(kp,cex=.8)
# kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
# for(i in names(genome_split)){
#    sub_cytoband <- cytoband_part_1[cytoband_part_1@seqnames==i]
#    sub_genome <- genome_part_1[genome_part_1@seqnames==i] 
#    sub_orientation = orientation_split[i][[1]]
#    sub_kprect = kprect_split[i][[1]]
#    plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
# }
# dev.off()


# genome_part_2 <- custom.genome[seqnames(custom.genome) %in% chrom_split_order$`1`]
# cytoband_part_2 <- subsetByOverlaps(custom.cytobands,genome_part_2)
# seqlevels(genome_part_2) <- seqlevelsInUse(genome_part_2)
# levels(genome_part_2) <- seqlevels(genome_part_2)
# seqlevels(cytoband_part_2) <- seqlevelsInUse(cytoband_part_2)
# levels(cytoband_part_2) <- seqlevels(cytoband_part_2)

# fig_out = paste(sample_handle, 'karyotype_split_2_of_2', '.pdf',sep='')
# pdf(fig_out, width=15,height=40,pointsize=1)
# genome_split = split(genome_part_2,genome_part_2@seqnames)
# kp <- plotKaryotype(chromosomes=chrom_split_order$`1`, genome=genome_part_2, plot.type=2, cytobands=cytoband_part_2, plot.params=pp, cex=3.5)
# kpAddBaseNumbers(kp,cex=.8)
# kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE)
# for(i in names(genome_split)){
#    sub_cytoband <- cytoband_part_2[cytoband_part_2@seqnames==i]
#    sub_genome <- genome_part_2[genome_part_2@seqnames==i] 
#    sub_orientation = orientation_split[i][[1]]
#    sub_kprect = kprect_split[i][[1]]
#    plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp)
# }
# dev.off()


x <- data.frame()
write.table(x, file=processed,col.names=FALSE)
