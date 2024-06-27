# Functions for plotting for Kazer et al., Immunity

require(dplyr)
require(ggplot2)

# Plots the frequency of each subpopultion over the total number of cells within a seurat meta.data table (or set by cell.denom). 
# Compared to function above, takes in a metadata table, adds zeros for missing cell clusters, no tissue is specified.
# input your column names without quotes

CalcFreqByHashFromMD = function(metadata, cell.denom = NULL, col = cluster.names, which.x.axis = timepoint, roi = NULL,
                                palette = NULL, plot.height = 10, plot.width = 10, plot.ncol = 3, filename = NULL){
  
  #if region is specified, subset the metadata table down to only that region
  if(!is.null(region)){
    metadata = subset(metadata, region == roi)
  }
  
  #remove any levels no longer remaining and any unassigned cells
  metadata$assignment = as.character(metadata$assignment)
  metadata = metadata[metadata$assignment != "",]
  metadata = metadata[metadata[[deparse(substitute(col))]] != "unknown", ]
  
  #create a list to fill later for plotting
  p.list = list()
  
  #Create frequency table for plotting
    # denominator is all cells from each replicate
  if(is.null(cell.denom)){ 
    freq.by.hash = metadata %>% dplyr::group_by({{which.x.axis}}, assignment, {{col}}) %>%
      dplyr::summarise(count = n()) %>% dplyr::mutate(freq = count/sum(count)) %>% ungroup() %>% 
      tidyr::complete(nesting(!!ensym(which.x.axis), assignment), {{col}}, fill = list(count = 0, freq = 0), explicit = FALSE) #this syntax is weird for nesting...
    
  } else { #choose the denominator (common)
    freq.by.hash = metadata %>% dplyr::group_by({{which.x.axis}}, assignment, {{col}}) %>%
      dplyr::summarise(count = n()) %>% dplyr::mutate(freq = count/cell.denom) %>% ungroup() %>% 
      tidyr::complete(nesting(!!ensym(which.x.axis), assignment), {{col}}, fill = list(count = 0, freq = 0), explicit = FALSE)
  }
  
  Cluster.Names = sort(unique(metadata[[deparse(substitute(col))]]))
  colnames(freq.by.hash)[3] = "Cluster.Names" 
    
  if(!is.null(filename)){
    #loop through the clusters for plotting
    for(i in 1:length(Cluster.Names)){
      plot.table = freq.by.hash[which(freq.by.hash$Cluster.Names == Cluster.Names[i]),]
      p.list[[i]] = ggplot(plot.table, aes(x={{which.x.axis}}, y=freq*100, color={{which.x.axis}})) + geom_point(size = 2.5) +
        scale_color_manual(values = palette) + 
        labs(x = "Time Point", y = "Frequency (%)", title = Cluster.Names[i]) + 
        guides(color = "none") + ylim(0,max(plot.table$freq)*100) +
        theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
    }
    plot_grid(plotlist = p.list, ncol = plot.ncol)
    ggsave(filename = filename, width = plot.width, height = plot.height)
  }
  
  return(freq.by.hash)
}


# Plots a dot plot over all timepoints with a smoothed loess line. If only one cluster is provided, points are colored by the given palette.
# If more than one cluster is provided, points are colored by the cluster they belong to with a matching colored line.
PlotTemportalAbundance = function(abundances, clusters, time.axis = timepoint, group = Cluster.Names, 
                                  value = log2CP3k, tp.palette = NULL, group.palette = NULL, 
                                  y.axis.title = "log-normalized abundance", span = 0.5){
  #browser()
  coi.table = abundances %>% dplyr::filter({{group}} %in% clusters)
  
  if(length(clusters) == 1){
    ggplot(coi.table, aes(x = {{time.axis}}, y = {{value}})) + geom_point(aes(col = {{time.axis}}), size=3) +
      geom_smooth(method = 'loess', mapping = aes(x = as.numeric({{time.axis}}), y = {{value}}), se = FALSE, span = span, col = "red") + 
      scale_colour_manual(values = tp.palette) + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) + 
      labs(x = "", y = y.axis.title, col = "Time Point", title = clusters) + ylim(0,NA)
  }
  else{
    ggplot(coi.table, aes(x = {{time.axis}}, y = {{value}}, col = {{group}})) + geom_point(aes(shape = {{group}}), size = 3) +
      geom_smooth(method = 'loess', mapping = aes(x = as.numeric({{time.axis}}), y = {{value}}), se = FALSE, span = span) + 
      scale_colour_manual(values = group.palette) + guides(shape = "none") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +  
      labs(x = "", y = y.axis.title, col = "Clusters", title = paste(clusters, collapse = " & ")) + ylim(0,NA)
  }
}


# Plots dot plot like above, but only for one cluster and separates the data into two (or more) time courses based on the `split.by` column
PlotTemporalAbundanceComparison = function(abundances, cluster, time.axis = timepoint.comp, group = Cluster.Names, 
                                           split.by = experiment, value = log2CP3k, tp.palette = NULL, group.palette = NULL, 
                                           y.axis.title = "log-normalized abundance", span = 0.5){
  
  coi.table = abundances %>% dplyr::filter({{group}} == cluster)
  
  ggplot(coi.table, aes(x = {{time.axis}}, y = {{value}}, col = {{split.by}})) + geom_jitter(aes(shape = {{split.by}}), size = 3, width = 0.1, height = 0) +
    geom_smooth(method = 'loess', mapping = aes(x = as.numeric({{time.axis}}), y = {{value}}), se = FALSE, span = span) + 
    scale_colour_manual(values = group.palette) + guides(shape = "none") + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +  
    labs(x = "", y = y.axis.title, col = "Experiment", title = cluster) + ylim(-0.5,NA)
}







# #Plots the frequency of each subpopultion over the total number of cells within srobj (or set by cell.denom). 
# #Specify tissue and x-axis for splitting of groups.
# 
# CalcFreqByHash = function(srobj, tissue = NULL, cell.denom = NULL, which.names = "Short", which.x.axis = "TimePoint",
#                           plot.height = 10, plot.width = 10, plot.ncol = 3, filename = NULL){
#   
#   #create a list to fill later for plotting
#   p.list = list()
#   
#   #short or detailed names
#   if(which.names == "Short"){
#     #Create frequency table for plotting
#     if(is.null(cell.denom)){ # denominator automatically calculated for each sample
#       freq.by.hash = srobj@meta.data %>% group_by(orig.ident, assignment, Short.Name) %>%
#         summarise(count = n(), .drop = FALSE) %>% mutate(freq = count/sum(count)) #%>% complete(assignment, fill = list(count = 0, freq = 0))
#     } else { #choose the denominator (common)
#       freq.by.hash = srobj@meta.data %>% group_by(orig.ident, assignment, Short.Name) %>%
#         summarise(count = n(), .drop = FALSE) %>% mutate(freq = count/cell.denom) #%>% complete(assignment, fill = list(count = 0, freq = 0))
#     }
#     Cluster.Names = levels(srobj$Short.Name)
#     colnames(freq.by.hash)[3] = "Cluster.Names" 
#     
#   } else if (which.names == "Detailed") {
#     if(is.null(cell.denom)){
#       freq.by.hash = srobj@meta.data %>% group_by(orig.ident, assignment, Detailed.Name, .drop = FALSE) %>%
#         summarise(count = n()) %>% mutate(freq = count/sum(count)) #%>% complete(assignment, fill = list(count = 0, freq = 0))
#     } else {
#       freq.by.hash = srobj@meta.data %>% group_by(orig.ident, assignment, Detailed.Name, .drop = FALSE) %>%
#         summarise(count = n()) %>% mutate(freq = count/cell.denom) #%>% complete(assignment, fill = list(count = 0, freq = 0))
#     }
#     Cluster.Names = levels(srobj$Detailed.Name)
#     colnames(freq.by.hash)[3] = "Cluster.Names" 
#   }  
#   
#   #loop through the clusters
#   for(i in 1:length(Cluster.Names)){
#     plot.table = freq.by.hash[intersect(which(freq.by.hash$Cluster.Names == Cluster.Names[i]), grep(tissue, freq.by.hash$assignment)),]
#     p.list[[i]] = ggplot(plot.table, aes(x=orig.ident, y=freq*100, color=orig.ident)) + geom_point(size = 2.5) +
#       scale_color_manual(values = sample.palette) + 
#       labs(x = "Time Point", y = "Frequency w/in Sample (%)", title = paste(Cluster.Names[i], tissue, sep = ", ")) + 
#       guides(color = "none") + ylim(0,max(plot.table$freq)*100) +
#       theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
#   }
#   plot_grid(plotlist = p.list, ncol = plot.ncol)
#   ggsave(filename = filename, width = plot.width, height = plot.height)
#   
#   return(freq.by.hash)
# }

