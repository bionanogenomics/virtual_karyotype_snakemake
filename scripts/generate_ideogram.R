args <- commandArgs(trailingOnly=TRUE)
print(args)
cyto <- args[1]
genome <- args[2]
orientation <- args[3]
kprect <- args[4]
processed <- args[5]
sample_handle <- args[6]
svlabel_handle<- args[7]

library(karyoploteR)

read_file_or_warn <- function(filepath) {
  if (file.info(filepath)$size == 1) {
    warning("The file is empty!")
    return(NULL)
  } else {
    return(read.table(filepath, sep='\t', header=T))
  }
}

custom.cytobands <- toGRanges(cyto)
custom.genome <- toGRanges(genome)
orientation_frame = read_file_or_warn(orientation)
kprect_frame = read.table(kprect, sep='\t', header=T)
svlabel_frame = read_file_or_warn(svlabel_handle)


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
   # Check if sub_orientation is NULL
   if (!is.null(sub_orientation)) {
      x0 = as.numeric(sub_orientation[2])
      x1 = as.numeric(sub_orientation[3])
      col = sub_orientation[20]
      kpArrows(kp, chr=path, x0=x0, x1=x1, y0=0.5, y1=0.5, data.panel=2, lwd=1, angle=20, col=col, length=0.5)
   }
}

plot_svlabel <- function(
   sub_svlabel, path, kp
){
   # Check if sub_orientation is NULL
   if (!is.null(sub_svlabel)) {
      print(sub_svlabel)
      labels = sub_svlabel[15]
      x = as.numeric(sub_svlabel[14])
      kpPlotMarkers(kp, chr=path, x=x, y=0.25, labels=labels, data.panel=2, cex=1.15)
   }
}

plot_ideogram <- function(
   path, fig_out, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper 
){
   pdf(fig_out)
   kp <- plotKaryotype(chromosomes=path, genome=sub_genome, plot.type=2, cytobands=sub_cytoband, plot.params=pp)
   kpAddBaseNumbers(kp)
   kpAddCytobandLabels(kp, force.all=TRUE, srt=90, col='#2C02FD', cex=0.5)
   apply(sub_kprect, 1, plot_kprects, path=path, color_mapper=color_mapper, kp=kp)
   if (!is.null(sub_orientation)){
      apply(sub_orientation, 1, plot_orientation, path=path, kp=kp)
   }
   dev.off()
}

plot_total_ideogram <- function(
   path, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp, sub_svlabel
){
   apply(sub_kprect, 1, plot_kprects, path=path, color_mapper=color_mapper, kp=kp)
   if (!is.null(sub_orientation)){
      apply(sub_orientation, 1, plot_orientation, path=path, kp=kp)
   }
   if (!is.null(sub_svlabel)){
      apply(sub_svlabel, 1, plot_svlabel, path=path, kp=kp)
   }
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

add_kp_labels <- function(karyoplot, chr.names=NULL, xoffset=0, yoffset=-30, ...) {
  #Validate parameters
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(is.null(chr.names)) chr.names <- karyoplot$chromosomes
  if(length(chr.names)==0) stop("In kpAddChromosomeNames: chr.names must have at least one element.")
  if(!all(methods::is(chr.names, "character"))) stop("In kpAddChromosomeNames: all elements of chr.names must be characters.")
  #Begin plotting
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  bb <- getChromosomeNamesBoundingBox(karyoplot)  
  x <- (bb$x0+bb$x1)/2 + xoffset
  y <- (bb$y0+bb$y1)/2 + yoffset
  for (name in names(x)) {
    if (grepl("dic|der|\\(t\\(", name)) {
      text(x = 0.29, y = y[[name]], label = "*", cex=20, srt=90, pos=4, col='red')
    }
  invisible(karyoplot)
  }
}



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


second_largest <- function(coord_frame) {
    # Sort the 'x' column in descending order
    sorted_values <- sort(coord_frame$x, decreasing = TRUE)
  
    # Return the second value from the sorted list
    return(sorted_values[2])
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
   max_x <- second_largest(coord_frame) + .035
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

calc_chrom_width <- function(pp) {
  # Below needs to be calculated for number of paths in coord frame 
  data1outmargin <- pp$data1outmargin 
  data2outmargin <- pp$data2outmargin 
  data1inmargin <- pp$data1inmargin
  data2inmargin <- pp$data2inmargin 
  ideogramheight <- pp$ideogramheight
  data1height <- pp$data1height
  data2height <- pp$data2height

  chrom_width <- data1outmargin + data2outmargin + data1inmargin + data2inmargin + ideogramheight + data1height + data2height
  return(chrom_width)
}

add_labels <- function(kp,cyto_first,pp) {
   coord_frame <- get_coords(kp,cyto_first)
   coord_frame$chrom <- droplevels(coord_frame$chrom)
   chrom_coords<-calculate_coordinate_space(coord_frame,pp)
   plot_labels_by_chrom(chrom_coords)   
}

split_by_column <- function(data_frame, column_name) {
  
  # Check if data frame is NULL
  if (is.null(data_frame)) {
    warning("The data frame is NULL!")
    return(NULL)
  }
  # Check if column exists in data frame
  if (!column_name %in% names(data_frame)) {
    warning(paste("Column", column_name, "not found in the data frame!"))
    return(NULL)
  }
  # Perform the split operation
  result <- split(data_frame, data_frame[[column_name]])
  return(result)
}

get_subframe <- function(frame_list, index) {
  # Check if the list is NULL or empty
  if (is.null(frame_list) || length(frame_list) == 0) {
    warning("The list is NULL or not populated!")
    return(NULL)
  }
  # Check if the index is valid
  if(!(i %in% names(frame_list))){
    warning(paste("Index", index, "is out of bounds!"))
    return(NULL)
  }
  # Return the subframe
  subframe <- frame_list[index][[1]]
  if (is.null(subframe) || nrow(subframe) == 0) {
    warning("The subframe is NULL or empty!")
    return(NULL)
  }
  return(subframe)
}

kpAddBaseNumbers <- function(karyoplot, tick.dist=20000000, tick.len=5, 
                             units="auto", add.units=FALSE,
                             digits=2, minor.ticks=TRUE, 
                            minor.tick.dist=5000000, minor.tick.len=2,  cex=0.5, 
                            tick.col=NULL, minor.tick.col=NULL, clipping=TRUE,  ...) {
  
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  if(!(units %in% c("auto", "b", "kb", "Kb", "mb", "Mb"))) stop("invalid units. Must be one of: 'auto', 'b', 'Kb' or 'Mb'")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
  
  
  toLabel <- function(n, units, add.units, digits) {
    if(add.units==TRUE) {
      unit.labels <- c("b", "Kb", "Mb")
    } else {
      unit.labels <- c("", "", "")
    }
    if(units == "auto") {
      if(abs(n)<1000) { units <- "b"}
      else if(abs(n)<1000000) {units <- "Kb"}
      else {units <- "Mb"}
    }
    if(tolower(units) == "b") return(paste0(as.character(n), unit.labels[1])) #b
    if(tolower(units) == "kb") return(paste0(as.character(round(n/1000, digits=digits)), unit.labels[2])) #Kb
    return(paste0(as.character(round(n/1000000, digits=digits)), unit.labels[3])) #Mb
  }
  
  old.scipen <- options("scipen")
  options(scipen=999)
  on.exit(options(scipen=old.scipen), add=TRUE)

 
  #For every chromsome
  chromosomes <- setdiff(kp$chromosomes,"SCALE")
  for(chr.name in chromosomes) {

    chr <- karyoplot$genome[chr.name]
    #Major ticks
    num.ticks <- width(chr)/tick.dist + 1
  
    tick.pos <- start(chr) + (tick.dist*(0:(num.ticks-1))) - 1 
    tick.pos[1] <- start(chr)
    
    #if zoomed in, keep only the ticks in the plot region
    if(karyoplot$zoom==TRUE) {
      if(clipping==TRUE) {
        tick.pos <- tick.pos[tick.pos >= start(karyoplot$plot.region) & tick.pos<= end(karyoplot$plot.region)]
      }
    }
    
    if(length(tick.pos)>0) {#We have to test here and cannot test on num.ticks to take the zooming into account
      tick.labels <- sapply(tick.pos, toLabel, units=units, add.units=add.units, digits=digits)
        #
      #function(x){return(toLabel(x, units=units, add.units=add.units, digits=digits))}
      
      
      xplot <- ccf(chr=rep(chr.name, length(tick.pos)), x=tick.pos, data.panel="ideogram")$x
      y0plot <- mids(chr.name)-karyoplot$plot.params$ideogramheight/2
      if(is.null(tick.col)) {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len, ...)
      } else {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len, col=tick.col, ...)
      }
      graphics::text(x=xplot, y=y0plot-tick.len, labels=tick.labels, pos=1, cex=cex, offset=0.1, ...)
    }
    #Minor ticks
    if(minor.ticks) {
      if(width(chr)>minor.tick.dist) {
        minor.num.ticks <- floor(width(chr)/minor.tick.dist)
        minor.tick.pos <- start(chr) + (minor.tick.dist*(seq_len(minor.num.ticks))) - 1
  
        #if zoomed in, keep only the ticks in the plot region
        if(karyoplot$zoom==TRUE) {
          if(clipping==TRUE) {
          minor.tick.pos <- minor.tick.pos[minor.tick.pos >= start(karyoplot$plot.region) & minor.tick.pos<= end(karyoplot$plot.region)]
          }
        }
        if(length(minor.tick.pos)>0) { 
          xplot <- ccf(chr=rep(chr.name, length(minor.tick.pos)), x=minor.tick.pos , data.panel="ideogram")$x
          y0plot <- mids(chr.name) - karyoplot$plot.params$ideogramheight/2
          if(is.null(minor.tick.col)) {
            graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len, ...)       
          } else {
            graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len, col=minor.tick.col, ...)       
          }
        }
      }
    }
  }
  
  invisible(karyoplot)
}


labelChromosomes <- function(kp, cex=1) {
  
  # Retrieve the r0 and r1 values for plotting area to position labels appropriately
  r0 <- min(kp$ideogram$ylim)
  r1 <- max(kp$ideogram$ylim)

  # A buffer to shift labels down so they aren't too close to the chromosomes
  y.buffer <- 0.02 * (r1 - r0)

  # Iterate through the chromosomes
  for(chrom in names(kp$chromosomes)) {
    
    # Check if there's a homologous pair
    if(sum(kp$ideogram$chromosome == chrom) == 2) {
      
      # If yes, find mid-point between homologous pair to place the label
      x.positions <- kp$ideogram$baseMid[kp$ideogram$chromosome == chrom]
      x.mid <- mean(x.positions)
      
    } else {
      
      # If no pair, place label directly below the chromosome
      x.mid <- kp$ideogram$baseMid[kp$ideogram$chromosome == chrom]
      
    }
    
    # Add the label to the plot
    kpText(kp, x=x.mid, y=r0 - y.buffer, labels=chrom, cex=cex, r0=r0, r1=r1)
  }
}

extract_acen <- function(df) {
  # Filter rows where gieStain equals 'acen'
  filtered_df <- df[df$gieStain == 'acen', c('seqnames', 'chrom')]
  
  # For each unique seqnames, select the first chrom
  selected_rows <- !duplicated(filtered_df$seqnames)
  unique_filtered_df <- filtered_df[selected_rows, ]
  
  # Setting row names
  rownames(unique_filtered_df) <- unique_filtered_df$seqnames
  
  return(unique_filtered_df)
}

add_sv_labels <- function(karyoplot, svlabel_frame, chr.names=NULL, xoffset=0, yoffset=-30, ...) {
  #Validate parameters
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  # If svlabel_frame is NULL, don't do anything further
  if(is.null(svlabel_frame)) {
    return(invisible(NULL))
  }

  if(is.null(chr.names)) chr.names <- karyoplot$chromosomes
  if(length(chr.names)==0) stop("In kpAddChromosomeNames: chr.names must have at least one element.")
  if(!all(methods::is(chr.names, "character"))) stop("In kpAddChromosomeNames: all elements of chr.names must be characters.")
  #Begin plotting
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  bb <- getChromosomeNamesBoundingBox(karyoplot)  
  x <- (bb$x0+bb$x1)/2 + xoffset
  y <- (bb$y0+bb$y1)/2 + yoffset
  for (name in names(x)) {
    if (name %in% svlabel_frame$chr) {
      text(x = 0.29, y = y[[name]], label = "*", cex=20, srt=90, pos=4, col='red')
    }
  invisible(karyoplot)
  }
}

pp <- getDefaultPlotParams(plot.type=2)
pp$leftmargin <- 0.3
pp$data2height <- 30

color_mapper <- mapColors(unique(custom.cytobands$chrom))
genome_split = split(custom.genome,custom.genome@seqnames)
cytoband_split = split(custom.cytobands, custom.cytobands@seqnames)
# orientation_split = split(orientation_frame, orientation_frame$chr)
orientation_split = split_by_column(orientation_frame, 'chr')
kprect_split = split(kprect_frame, kprect_frame$chr)
svlabel_split = split_by_column(svlabel_frame, 'chr')

# for(i in names(genome_split)){
#    path_handle = strsplit(i, " ")[[1]][1]
#    fig_out = paste(sample_handle, path_handle, '.pdf',sep='')
#    sub_cytoband <- custom.cytobands[custom.cytobands@seqnames==i]
#    sub_genome <- custom.genome[custom.genome@seqnames==i] 
#    # sub_orientation = orientation_split[i][[1]]
#    sub_orientation = get_subframe(orientation_split, i)
#    sub_kprect = kprect_split[i][[1]]
#    plot_ideogram(i, fig_out, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper)
# }

################
pp <- getDefaultPlotParams(plot.type=2)
pp$leftmargin <- 0.3
pp$data2height <- 100
pp$ideogramheight <- 80 
pp$dataideogrammin <- -1
pp$dataideogrammax <- 1
pp$data2inmargin <- 20
pp$data1inmargin <- 20
pp$data1height <- 100
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
custom.genome@seqinfo@seqlengths <- width(custom.genome)
seqinfo_obj <- seqinfo(custom.genome)
max_seqname <- seqnames(seqinfo_obj)[which.max(seqlengths(seqinfo_obj))]

subset_gr_max <- custom.genome[seqnames(custom.genome) == max_seqname]
seqlevels(subset_gr_max) <- seqlevelsInUse(subset_gr_max)
gr_scaler <- GRanges(seqnames = c('SCALE'),ranges=subset_gr_max@ranges[1])


contig_split <- round(length(chrom_order)/2)

chrom_split_order <- chunk.2(chrom_order, 2, force.number.of.groups=T)

genome_part_1 <- custom.genome[seqnames(custom.genome) %in% chrom_split_order$`2`]
cytoband_part_1 <- subsetByOverlaps(custom.cytobands,genome_part_1)
seqlevels(genome_part_1) <- seqlevelsInUse(genome_part_1)
levels(genome_part_1) <- seqlevels(genome_part_1)
seqlevels(cytoband_part_1) <- seqlevelsInUse(cytoband_part_1)
levels(cytoband_part_1) <- seqlevels(cytoband_part_1)


cyto_first = extract_acen(cyto)
#cyto_first <- get_first_rows_per_seqnames(sub_cyto)

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


kpAddChromosomeNames <- function(karyoplot, chr.names=NULL, xoffset=0, yoffset=0, ...) {
  
  # Validate parameters
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(is.null(chr.names)) {
    chr.names <- karyoplot$chromosomes
  }

  # Ensure "SCALE" is excluded
  chr.names <- setdiff(chr.names, "SCALE")

  if(length(chr.names)==0) stop("In kpAddChromosomeNames: chr.names must have at least one element.")
  if(!all(methods::is(chr.names, "character"))) stop("In kpAddChromosomeNames: all elements of chr.names must be characters.")

  # Begin plotting
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  bb <- getChromosomeNamesBoundingBox(karyoplot)
  
  x <- (bb$x0+bb$x1)/2 + xoffset
  y <- (bb$y0+bb$y1)/2 + yoffset
  
  # Filter out the positions related to "SCALE" chromosome
  x <- x[chr.names]
  y <- y[chr.names]
  
  graphics::text(x=x, y=y, labels=chr.names, ...)
  
  invisible(karyoplot)
}


fig_out = paste(sample_handle, 'karyotype_split_1_of_4', '.pdf',sep='')
pdf(fig_out, width=20, height=40, pointsize=1)
genome_split = split(genome_part1_5, genome_part1_5@seqnames)
# Exclude the 'SCALE' chromosome using the zoom parameter
visible_chromosomes <- setdiff(levels(genome_part1_5@seqnames), "SCALE")
kp <- plotKaryotype(chromosomes = c(visible_chromosomes, 'SCALE'),
                    genome = c(genome_part1_5, gr_scaler),
                    plot.type = 2,
                    cytobands = cytoband_part1_5,
                    plot.params = pp,
                    lwd = 0.1,
                    cex = 3.5,
                    labels.plotter=kpAddChromosomeNames)

# Ensure other plotting functions only act on visible chromosomes
for(i in visible_chromosomes) {
   kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", chromosomes = i)
   kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE, chromosomes = i)
   sub_cytoband <- cytoband_part1_5[cytoband_part1_5@seqnames == i]
   sub_genome <- genome_part1_5[genome_part1_5@seqnames == i] 
   sub_orientation = get_subframe(orientation_split, i)
   sub_kprect = kprect_split[i][[1]]
   sub_svlabel = get_subframe(svlabel_split, i)
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp, sub_svlabel)
}
add_labels(kp, cyto_first, pp)
add_kp_labels(kp)
add_sv_labels(kp, svlabel_frame)
dev.off()

fig_out = paste(sample_handle, 'karyotype_split_2_of_4', '.pdf',sep='')
pdf(fig_out, width=20, height=40, pointsize=1)
genome_split = split(genome_part6_12, genome_part6_12@seqnames)
# Exclude the 'SCALE' chromosome using the zoom parameter
visible_chromosomes <- setdiff(levels(genome_part6_12@seqnames), "SCALE")
kp <- plotKaryotype(chromosomes = c(visible_chromosomes, 'SCALE'),
                    genome = c(genome_part6_12, gr_scaler),
                    plot.type = 2,
                    cytobands = cytoband_part6_12,
                    plot.params = pp,
                    lwd = 0.1,
                    cex = 3.5,
                    labels.plotter=kpAddChromosomeNames)

# Ensure other plotting functions only act on visible chromosomes
for(i in visible_chromosomes) {
   kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", chromosomes = i)
   kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE, chromosomes = i)
   sub_cytoband <- cytoband_part6_12[cytoband_part6_12@seqnames == i]
   sub_genome <- genome_part6_12[genome_part6_12@seqnames == i] 
   sub_orientation = get_subframe(orientation_split, i)
   sub_kprect = kprect_split[i][[1]]
   sub_svlabel = get_subframe(svlabel_split, i)
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp, sub_svlabel)
}
add_labels(kp, cyto_first, pp)
add_kp_labels(kp)
add_sv_labels(kp, svlabel_frame)
dev.off()

fig_out = paste(sample_handle, 'karyotype_split_3_of_4', '.pdf',sep='')
pdf(fig_out, width=20, height=40, pointsize=1)
genome_split = split(genome_part13_18, genome_part13_18@seqnames)
# Exclude the 'SCALE' chromosome using the zoom parameter
visible_chromosomes <- setdiff(levels(genome_part13_18@seqnames), "SCALE")
kp <- plotKaryotype(chromosomes = c(visible_chromosomes, 'SCALE'),
                    genome = c(genome_part13_18, gr_scaler),
                    plot.type = 2,
                    cytobands = cytoband_part13_18,
                    plot.params = pp,
                    lwd = 0.1,
                    cex = 3.5,
                    labels.plotter=kpAddChromosomeNames)

# Ensure other plotting functions only act on visible chromosomes
for(i in visible_chromosomes) {
   kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", chromosomes = i)
   kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE, chromosomes = i)
   sub_cytoband <- cytoband_part13_18[cytoband_part13_18@seqnames == i]
   sub_genome <- genome_part13_18[genome_part13_18@seqnames == i] 
   sub_orientation = get_subframe(orientation_split, i)
   sub_kprect = kprect_split[i][[1]]
   sub_svlabel = get_subframe(svlabel_split, i)
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp, sub_svlabel)
}
add_labels(kp, cyto_first, pp)
add_kp_labels(kp)
add_sv_labels(kp, svlabel_frame)
dev.off()


fig_out = paste(sample_handle, 'karyotype_split_4_of_4', '.pdf',sep='')
pdf(fig_out, width=20, height=40, pointsize=1)
genome_split = split(genome_part19_Y, genome_part19_Y@seqnames)
# Exclude the 'SCALE' chromosome using the zoom parameter
visible_chromosomes <- setdiff(levels(genome_part19_Y@seqnames), "SCALE")
kp <- plotKaryotype(chromosomes = c(visible_chromosomes, 'SCALE'),
                    genome = c(genome_part19_Y, gr_scaler),
                    plot.type = 2,
                    cytobands = cytoband_part19_Y,
                    plot.params = pp,
                    lwd = 0.1,
                    cex = 3.5,
                    labels.plotter=kpAddChromosomeNames)

# Ensure other plotting functions only act on visible chromosomes
for(i in visible_chromosomes) {
   kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray", chromosomes = i)
   kpAddCytobandLabels(kp, srt=90, col='#2C02FD', cex=1.5, force.all=TRUE, chromosomes = i)
   sub_cytoband <- cytoband_part19_Y[cytoband_part19_Y@seqnames == i]
   sub_genome <- genome_part19_Y[genome_part19_Y@seqnames == i] 
   sub_orientation = get_subframe(orientation_split, i)
   sub_kprect = kprect_split[i][[1]]
   sub_svlabel = get_subframe(svlabel_split, i)
   plot_total_ideogram(i, sub_genome, sub_cytoband, sub_orientation, sub_kprect, color_mapper, kp, sub_svlabel)
}
add_labels(kp, cyto_first, pp)
add_kp_labels(kp)
add_sv_labels(kp, svlabel_frame)
dev.off()


x <- data.frame()
write.table(x, file=processed,col.names=FALSE)
