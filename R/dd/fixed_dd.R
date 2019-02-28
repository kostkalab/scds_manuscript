
#- doubletDecon: Synthetic_Doublets.R @ line 76 uses a parameter 'location'
#                that is not defined in this function and not a parameter.
#                The only change in the code below is to  make it an explicit
#                parameter
#
#                VERSION DD: Latest commit 606ec17  on Jan 9
#                assessed on: 2/24/2019

mySynthetic_Doublets <- function (data, groups, groupsMedoids, newMedoids, num_doubs,
          log_file_name, only50, location)
{
  groups = groupsMedoids
  ndub = num_doubs
  pairs = combn(unique(groups[, 2]), 2)
  doubletCellsInput = as.data.frame(matrix(ncol = 2, nrow = ((length(pairs)/2) *
                                                               ndub)))
  doubletCellsInput2 = as.data.frame(matrix(ncol = 4, nrow = ((length(pairs)/2) *
                                                                ndub)))
  for (pair in 0:((length(pairs)/2) - 1)) {
    for (synth in 1:ndub) {
      doubletCellsInput[synth + (ndub * pair), 1] = sample(row.names(subset(groups,
                                                                            groups[, 2] == pairs[1, (pair + 1)])), 1, replace = FALSE)
      doubletCellsInput2[synth + (num_doubs * pair), 1] = doubletCellsInput[synth +
                                                                              (num_doubs * pair), 1]
      doubletCellsInput[synth + (ndub * pair), 2] = sample(row.names(subset(groups,
                                                                            groups[, 2] == pairs[2, (pair + 1)])), 1, replace = FALSE)
      doubletCellsInput2[synth + (num_doubs * pair), 2] = doubletCellsInput[synth +
                                                                              (num_doubs * pair), 2]
      doubletCellsInput2[synth + (num_doubs * pair), 3] = pairs[1,
                                                                (pair + 1)]
      doubletCellsInput2[synth + (num_doubs * pair), 4] = pairs[2,
                                                                (pair + 1)]
    }
  }
  doubletAverages = as.data.frame(matrix(ncol = nrow(doubletCellsInput),
                                         nrow = nrow(data)))
  row.names(doubletAverages) = row.names(data)
  doubletAverages[1, ] = rep((length(unique(groups[, 1])) +
                                1), ncol(doubletAverages))
  doubletAverages = as.data.frame(matrix(ncol = nrow(doubletCellsInput) *
                                           3, nrow = nrow(data)))
  row.names(doubletAverages) = row.names(data)
  doubletAverages[1, ] = rep((length(unique(groups[, 1])) +
                                1), ncol(doubletAverages) * 3)
  for (doublet in 1:(ncol(doubletAverages)/3)) {
    cell1 = as.character(doubletCellsInput[doublet, 1])
    cell2 = as.character(doubletCellsInput[doublet, 2])
    expression1 = data[2:nrow(data), which(colnames(data) ==
                                             cell1)]
    expression2 = data[2:nrow(data), which(colnames(data) ==
                                             cell2)]
    temp = cbind(expression1, expression2)
    newExpression = apply(temp, 1, weighted.mean, c(0.5,
                                                    0.5))
    newExpression_a = apply(temp, 1, weighted.mean, c(0.7,
                                                      0.3))
    newExpression_b = apply(temp, 1, weighted.mean, c(0.3,
                                                      0.7))
    doubletAverages[2:nrow(doubletAverages), doublet] = newExpression
    doubletAverages[2:nrow(doubletAverages), doublet + (ncol(doubletAverages)/3)] = newExpression_a
    doubletAverages[2:nrow(doubletAverages), doublet + ((ncol(doubletAverages)/3) *
                                                          2)] = newExpression_b
    colnames(doubletAverages)[doublet] = paste0(cell1, "-",
                                                cell2, "-even")
    colnames(doubletAverages)[doublet + (ncol(doubletAverages)/3)] = paste0(cell1,
                                                                            "-", cell2, "-one")
    colnames(doubletAverages)[doublet + ((ncol(doubletAverages)/3) *
                                           2)] = paste0(cell1, "-", cell2, "-two")
  }
  if (only50 == TRUE) {
    mult = 1
  }
  else {
    mult = 3
  }
  results = DeconRNASeq(doubletAverages[2:nrow(doubletAverages),
                                        ], newMedoids)
  resultsreadable = round(results$out.all * 100, 2)
  write.table(resultsreadable, paste0(location, "resultsreadable_synths.txt"),
              sep = "\t")
  row.names(resultsreadable) = colnames(doubletAverages)
  averagesAverages = as.data.frame(matrix(ncol = ncol(resultsreadable),
                                          nrow = (length(pairs)/2) * mult))
  colnames(averagesAverages) = colnames(resultsreadable)
  i = 1
  for (clust in 1:nrow(averagesAverages)) {
    averagesAverages[clust, ] = apply(resultsreadable[i:(i +
                                                           (num_doubs - 1)), ], 2, mean)
    if (clust %in% 1:(length(pairs)/2)) {
      row.names(averagesAverages)[clust] = paste0(pairs[1,
                                                        clust], "-", pairs[2, clust], "-even")
    }
    else if (clust %in% ((length(pairs)/2) + 1):((length(pairs)/2) *
                                                 2)) {
      row.names(averagesAverages)[clust] = paste0(pairs[1,
                                                        clust - (length(pairs)/2)], "-", pairs[2, clust -
                                                                                                 (length(pairs)/2)], "-one")
    }
    else {
      row.names(averagesAverages)[clust] = paste0(pairs[1,
                                                        clust - length(pairs)], "-", pairs[2, clust -
                                                                                             length(pairs)], "-two")
    }
    i = i + num_doubs
  }
  return(list(averagesAverages = averagesAverages, doubletCellsInput2 = doubletCellsInput2))
}

myMain_Doublet_Decon <- function (rawDataFile, groupsFile, filename, location, fullDataFile = NULL,
          removeCC = FALSE, species = "mmu", rhop = 1, write = TRUE,
          PMF = TRUE, useFull = FALSE, heatmap = TRUE, centroids = FALSE,
          num_doubs = 30, downsample = "none", sample_num = NULL, only50 = TRUE,
          min_uniq = 4)
{
  require(DeconRNASeq)
  require(gplots)
  require(dplyr)
  require(MCL)
  require(clusterProfiler)
  require(mygene)
  log_file_name = paste0(location, Sys.time(), ".log")
  log_con <- file(log_file_name)
  cat(paste0("filename: ", filename), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("location: ", location), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("removeCC: ", removeCC), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("species: ", species), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("rhop: ", rhop), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("write: ", write), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("PMF: ", PMF), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("useFull: ", useFull), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("heatmap: ", heatmap), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("centroids: ", centroids), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("num_doubs: ", num_doubs), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("downsample: ", downsample), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("sample_num: ", sample_num), file = log_file_name,
      append = TRUE, sep = "\n")
  cat(paste0("only50: ", only50), file = log_file_name, append = TRUE,
      sep = "\n")
  cat(paste0("min_uniq: ", min_uniq), file = log_file_name,
      append = TRUE, sep = "\n")
  if (is.character(rawDataFile) != TRUE & is.data.frame(rawDataFile) !=
      TRUE) {
    print("ERROR: rawDataFile must be a character string!")
  }
  if (is.character(groupsFile) != TRUE & is.data.frame(groupsFile) !=
      TRUE) {
    print("ERROR: groupsFile must be a character string!")
  }
  if (is.character(filename) != TRUE) {
    print("ERROR: filename must be a character string!")
  }
  if (is.character(location) != TRUE) {
    print("ERROR: location must be a character string!")
  }
  if (is.character(fullDataFile) != TRUE & is.null(fullDataFile) !=
      TRUE & is.data.frame(fullDataFile) != TRUE) {
    print("ERROR: fullDataFile must be a character string or NULL!")
  }
  if (is.logical(removeCC) != TRUE) {
    print("ERROR: removeCC must be TRUE or FALSE!")
  }
  if (is.character(species) != TRUE) {
    print("ERROR: species must be a character string!")
  }
  if (is.numeric(rhop) != TRUE) {
    print("ERROR: rhop must be numeric!")
  }
  if (is.logical(write) != TRUE) {
    print("ERROR: write must be TRUE or FALSE!")
  }
  if (is.logical(PMF) != TRUE) {
    print("ERROR: PMF must be TRUE or FALSE!")
  }
  if (is.logical(useFull) != TRUE) {
    print("ERROR: useFull must be TRUE or FALSE!")
  }
  if (is.logical(heatmap) != TRUE) {
    print("ERROR: heatmap must be TRUE or FALSE!")
  }
  if (is.logical(centroids) != TRUE) {
    print("ERROR: centroids must be TRUE or FALSE!")
  }
  if (is.numeric(num_doubs) != TRUE) {
    print("ERROR: numdoubs must be numeric!")
  }
  if (is.character(downsample) != TRUE) {
    print("ERROR: downsample must be a character string!")
  }
  if (is.numeric(sample_num) != TRUE & is.null(sample_num) !=
      TRUE) {
    print("ERROR: sample_num must be numeric or NULL!")
  }
  if (is.logical(only50) != TRUE) {
    print("ERROR: only50 must be TRUE or FALSE!")
  }
  if (is.numeric(min_uniq) != TRUE) {
    print("ERROR: min_uniq must be numeric!")
  }
  cat("Reading data...", file = log_file_name, append = TRUE,
      sep = "\n")
  cat("Reading data...", sep = "\n")
  if (class(rawDataFile) == "character") {
    rawData = read.table(rawDataFile, sep = "\t", header = T,
                         row.names = 1)
  }
  else {
    rawData = rawDataFile
  }
  if (class(groupsFile) == "character") {
    groups = read.table(groupsFile, sep = "\t", header = F,
                        row.names = 1)
  }
  else {
    groups = groupsFile
  }
  cat("Processing raw data...", file = log_file_name, append = TRUE,
      sep = "\n")
  cat("Processing raw data...", sep = "\n")
  data = Clean_Up_Input(rawData, groups, log_file_name = log_file_name)
  og_processed_data = data$processed
  groups = data$groups
  if (centroids == TRUE) {
    centroid_flag = TRUE
  }
  else {
    centroid_flag = FALSE
  }
  if (heatmap == TRUE) {
    cat("Creating original data heatmap...", file = log_file_name,
        append = TRUE, sep = "\n")
    cat("Creating original data heatmap...", sep = "\n")
    breaks = seq(0, as.numeric(quantile(data.matrix(data$processed[2:nrow(data$processed),
                                                                   2:ncol(data$processed)]), 0.99)), by = 0.05)
    mycol <- colorpanel(n = length(breaks) - 1, low = "black",
                        high = "yellow")
    suppressWarnings(DDheatmap(data.matrix(data$processed[2:nrow(data$processed),
                                                          2:ncol(data$processed)]), Colv = FALSE, Rowv = FALSE,
                               dendrogram = "none", col = mycol, ColSideColors = as.color(Renumber(data$processed[1,
                                                                                                                  2:ncol(data$processed)]), alpha = 1, seed = 4),
                               RowSideColors = as.color(Renumber(data$processed[2:nrow(data$processed),
                                                                                1]), alpha = 1, seed = 2), breaks = breaks, trace = "none",
                               na.rm = TRUE, margins = c(5, 5), labRow = NA, labCol = NA,
                               xlab = "Samples", ylab = "Genes", main = paste0("Original data: ",
                                                                               filename)))
  }
  if (removeCC == TRUE) {
    cat("Removing cell cycle clusters...", file = log_file_name,
        append = TRUE, sep = "\n")
    cat("Removing cell cycle clusters...", sep = "\n")
    data = Remove_Cell_Cycle(data$processed, species, log_file_name)
  }
  else {
    data = data$processed
  }
  if (write == TRUE) {
    write.table(data, paste0(location, "data_processed_",
                             filename, ".txt"), sep = "\t")
    write.table(groups, paste0(location, "groups_processed_",
                               filename, ".txt"), sep = "\t")
  }
  cat("Combining similar clusters...", file = log_file_name,
      append = TRUE, sep = "\n")
  cat("Combining similar clusters...", sep = "\n")
  BL = Blacklist_Groups(data, groups, rhop, centroid_flag,
                        log_file_name)
  newMedoids = BL$newMedoids
  groupsMedoids = BL$newGroups
  cat("Creating synthetic doublet profiles...", file = log_file_name,
      append = TRUE, sep = "\n")
  cat("Creating synthetic doublet profiles...", sep = "\n")
  if (.Platform$OS.type == "unix") {
    sink("/dev/null")
    synthProfilesx = mySynthetic_Doublets(data, groups, groupsMedoids,
                                        newMedoids, num_doubs, log_file_name = log_file_name,
                                        only50 = only50,
                                        location)
    sink()
  }
  else {
    synthProfilesx = mySynthetic_Doublets(data, groups, groupsMedoids,
                                        newMedoids, num_doubs, log_file_name = log_file_name,
                                        only50 = only50,
                                        location)
  }
  synthProfiles = synthProfilesx$averagesAverages
  doubletCellsInput2 = synthProfilesx$doubletCellsInput2
  if (write == TRUE) {
    write.table(doubletCellsInput2, paste0(location, "Synth_doublet_info_",
                                           filename, ".txt"), sep = "\t")
  }
  cat("Step 1: Removing possible doublets...", file = log_file_name,
      append = TRUE, sep = "\n")
  cat("Step 1: Removing possible doublets...", sep = "\n")
  if (.Platform$OS.type == "unix") {
    sink("/dev/null")
    doubletTable = Is_A_Doublet(data, newMedoids, groups,
                                synthProfiles, log_file_name = log_file_name)
    sink()
  }
  else {
    doubletTable = Is_A_Doublet(data, newMedoids, groups,
                                synthProfiles, log_file_name = log_file_name)
  }
  if (write == TRUE) {
    write.table(doubletTable$isADoublet, paste0(location,
                                                "DRS_doublet_table_", filename, ".txt"), sep = "\t")
    write.table(doubletTable$resultsreadable, paste0(location,
                                                     "DRS_results_", filename, ".txt"), sep = "\t")
  }
  cat("Step 2: Re-clustering possible doublets...", file = log_file_name,
      append = TRUE, sep = "\n")
  cat("Step 2: Re-clustering possible doublets...", sep = "\n")
  reclusteredData = Recluster(isADoublet = doubletTable$isADoublet,
                              data, groups, log_file_name = log_file_name)
  data = reclusteredData$newData2$processed
  groups = reclusteredData$newData2$groups
  if (write == TRUE) {
    write.table(data, paste0(location, "data_processed_reclust_",
                             filename, ".txt"), sep = "\t")
    write.table(groups, paste0(location, "groups_processed_reclust_",
                               filename, ".txt"), sep = "\t")
  }
  if (PMF == FALSE) {
    cat("SKIPPING Step 3: Rescuing cells with unique gene expression...",
        file = log_file_name, append = TRUE, sep = "\n")
    cat("SKIPPING Step 3: Rescuing cells with unique gene expression...",
        sep = "\n")
    PMFresults = NULL
  }
  else {
    cat("Step 3: Rescuing cells with unique gene expression...",
        file = log_file_name, append = TRUE, sep = "\n")
    cat("Step 3: Rescuing cells with unique gene expression...",
        sep = "\n")
    if (useFull == TRUE) {
      if (class(fullDataFile) == "character") {
        full_data = read.table(fullDataFile, sep = "\t",
                               header = T, row.names = 1)
      }
      else {
        full_data = fullDataFile
      }
      full_data2 = Clean_Up_Input(full_data, groups)$processed
      PMFresults = Pseudo_Marker_Finder(groups, data, full_data2,
                                        downsample, sample_num, min_uniq = min_uniq,
                                        log_file_name = log_file_name)
    }
    else {
      PMFresults = Pseudo_Marker_Finder(groups, data, full_data2 = NULL,
                                        downsample, sample_num, min_uniq = min_uniq,
                                        log_file_name = log_file_name)
    }
    if (write == TRUE) {
      write.table(PMFresults, paste0(location, "new_PMF_results_",
                                     filename, ".txt"), sep = "\t")
    }
  }
  allClusters = unique(groups[, 1])
  if (PMF == FALSE) {
    newDoubletClusters = allClusters
  }
  else {
    hallmarkClusters = as.numeric(unique(PMFresults[, 2]))
    newDoubletClusters = setdiff(allClusters, hallmarkClusters)
  }
  uniqueClusters = as.character(unique(groups[, 2]))
  DeconCalledFreq = as.data.frame(matrix(nrow = length(allClusters),
                                         ncol = 1), row.names = uniqueClusters)
  for (clus in 1:length(allClusters)) {
    temp1 = subset(doubletTable$isADoublet, Group_Cluster ==
                     uniqueClusters[clus])
    if (nrow(temp1) == 0) {
      DeconCalledFreq[clus, 1] = 100
    }
    else {
      DeconCalledFreq[clus, 1] = (length(which(temp1$isADoublet ==
                                                 TRUE))/nrow(temp1)) * 100
    }
  }
  if (PMF == FALSE) {
    finalDoubletClusters = which(DeconCalledFreq > 50)
  }
  else {
    finalDoubletClusters = intersect(which(DeconCalledFreq >
                                             50), newDoubletClusters)
  }
  finalDoubletCellCall = subset(groups, groups[, 1] %in% finalDoubletClusters)
  finalNotDoubletCellCall = subset(groups, !(groups[, 1] %in%
                                               finalDoubletClusters))
  if (write == TRUE) {
    write.table(finalDoubletCellCall, paste0(location, "Final_doublets_groups_",
                                             filename, ".txt"), sep = "\t")
    write.table(finalNotDoubletCellCall, paste0(location,
                                                "Final_nondoublets_groups_", filename, ".txt"), sep = "\t")
  }
  doublets_matrix = cbind(og_processed_data[, 1], og_processed_data[,
                                                                    which(colnames(og_processed_data) %in% row.names(finalDoubletCellCall))])
  if (write == TRUE) {
    write.table(doublets_matrix, paste0(location, "Final_doublets_exp_",
                                        filename, ".txt"), sep = "\t")
  }
  if (heatmap == TRUE) {
    cat("Creating doublets heatmap...", file = log_file_name,
        append = TRUE, sep = "\n")
    cat("Creating doublets heatmap...", sep = "\n")
    breaks = seq(0, as.numeric(quantile(data.matrix(doublets_matrix[2:nrow(doublets_matrix),
                                                                    2:ncol(doublets_matrix)]), 0.99)), by = 0.05)
    mycol <- colorpanel(n = length(breaks) - 1, low = "black",
                        high = "yellow")
    suppressWarnings(DDheatmap(data.matrix(doublets_matrix[2:nrow(doublets_matrix),
                                                           2:ncol(doublets_matrix)]), Colv = FALSE, Rowv = FALSE,
                               col = mycol, dendrogram = "none", ColSideColors = as.color(Renumber(doublets_matrix[1,
                                                                                                                   2:ncol(doublets_matrix)]), alpha = 1, seed = 4),
                               RowSideColors = as.color(Renumber(doublets_matrix[2:nrow(doublets_matrix),
                                                                                 1]), alpha = 1, seed = 2), breaks = breaks, trace = "none",
                               na.rm = TRUE, margins = c(5, 5), labRow = NA, labCol = NA,
                               xlab = "Samples", ylab = "Genes", main = paste0("Doublets: ",
                                                                               filename)))
  }
  nondoublets_matrix = cbind(og_processed_data[, 1], og_processed_data[,
                                                                       which(colnames(og_processed_data) %in% row.names(finalNotDoubletCellCall))])
  if (write == TRUE) {
    write.table(nondoublets_matrix, paste0(location, "Final_nondoublets_exp_",
                                           filename, ".txt"), sep = "\t")
  }
  if (heatmap == TRUE) {
    cat("Creating non-doublets heatmap...", file = log_file_name,
        append = TRUE, sep = "\n")
    cat("Creating non-doublets heatmap...", sep = "\n")
    breaks = seq(0, as.numeric(quantile(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix),
                                                                       2:ncol(nondoublets_matrix)]), 0.99)), by = 0.05)
    mycol <- colorpanel(n = length(breaks) - 1, low = "black",
                        high = "yellow")
    suppressWarnings(DDheatmap(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix),
                                                              2:ncol(nondoublets_matrix)]), Colv = FALSE, Rowv = FALSE,
                               col = mycol, dendrogram = "none", ColSideColors = as.color(Renumber(nondoublets_matrix[1,
                                                                                                                      2:ncol(nondoublets_matrix)]), alpha = 1, seed = 4),
                               RowSideColors = as.color(Renumber(nondoublets_matrix[2:nrow(nondoublets_matrix),
                                                                                    1]), alpha = 1, seed = 2), breaks = breaks, trace = "none",
                               na.rm = TRUE, margins = c(5, 5), labRow = NA, labCol = NA,
                               xlab = "Samples", ylab = "Genes", main = paste0("Non-Doublets: ",
                                                                               filename)))
  }
  close(log_con)
  return(list(data_processed = data, groups_processed = groups,
              DRS_doublet_table = doubletTable$isADoublet, DRS_results = doubletTable$resultsreadable,
              PMF_results = PMFresults, Final_doublets_groups = finalDoubletCellCall,
              Final_nondoublets_groups = finalNotDoubletCellCall, Final_doublets_exp = doublets_matrix,
              Final_nondoublets_exp = nondoublets_matrix, Synth_doublet_info = doubletCellsInput2))
}
