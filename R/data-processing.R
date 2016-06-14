## # Download the latest version of the timetree of life
getTTOL <- function(treeNumber=2, list.only=FALSE){
  timetrees <- c(
    "TTOL_all_unsmoothed",
    "TTOL_all_smoothed",
    "TTOL_animals_unsmoothed",
    "TTOL_plants_unsmoothed",
    "TTOL_fungi_unsmoothed",
    "TTOL_birds_unsmoothed",
    "TTOL_birds_smoothed_interpolated",
    "TTOL_mammals_unsmoothed",
    "TTOL_mammals_smoothed_interpolated",
    "TTOL_squamates_unsmoothed",
    "TTOL_amphibians_unsmoothed",
    "TTOL_actinopterygians_unsmoothed",
    "TTOL_arthropods_unsmoothed",
    "TTOL_mollusks_unsmoothed",
    "TTOL_cnidarians_unsmoothed",
    "TTOL_ascomycetes_unsmoothed",
    "TTOL_basidiomycetes_unsmoothed",
    "TTOL_angiosperms_unsmoothed",
    "TTOL_gymnosperms_unsmoothed",
    "TTOL_moniliformopses_unsmoothed",
    "TTOL_liverworts_unsmoothed",
    "TTOL_mosses_unsmoothed",
    "TTOL_chlorophytes_unsmoothed",
    "TTOL_diatoms_unsmoothed",
    "TTOL_phaeophyceans_unsmoothed",
    "TTOL_eubacteria_unsmoothed",
    "TTOL_archaebacteria_unsmoothed")
  if(list.only){
    print(paste(1:length(timetrees), timetrees, sep="."))
  } else {
    treename <- paste(treeNumber, timetrees[treeNumber], sep=".")
    url <- paste("http://www.biodiversitycenter.org/files/ttol/",treename,".nwk", sep="")
    tree <- read.tree(url)
    return(tree)
  }
}

## # Build taxonomy for the timetree of life
getOttIds <- function(taxalist, ncores=1, context=NULL){
  scipen <- options()$scipen
  digits <- options()$digits
  options("scipen"=100, "digits"=4)
  .taxalist <- gsub("_", " ", taxalist)
  .taxalist <- gsub(" sp.", "", .taxalist)
  tax <- parallel::mclapply(1:length(taxalist),  function(i) try(rotl::tnrs_match_names(.taxalist[i], do_approximate_matching =FALSE, context_name = context)), mc.cores=ncores)
  failed <- which(sapply(tax,function(x) class(x)[1]=="try-error"))
  if(length(failed)>0){
    tax[failed] <- parallel::mclapply(failed,  function(i) try(rotl::tnrs_match_names(.taxalist[i], do_approximate_matching =TRUE, context_name = context)), mc.cores=ncores)
  }
  stillfailed <- which(sapply(tax,function(x) if(class(x)[1]=="try-error"){TRUE} else {is.na(x$ott_id)} ))
  if(length(stillfailed>0)){
    tax[stillfailed] <- lapply(stillfailed, function(x) data.frame(search_string=.taxalist[x], unique_name=.taxalist[x], approximate_match=NA, ott_id=NA, is_synonym=NA, flags=NA, number_matches=0))
  }
  tax <- do.call(rbind, tax)
  genspec <- unname(sapply(tax[,2], function(x) paste(strsplit(x, split=" ")[[1]][1:2],collapse=" ")))
  genspec <- gsub(" (genus", " sp.", genspec, fixed=TRUE)
  genspec <- gsub(" NA", " sp.", genspec, fixed=TRUE)
  if(sum(duplicated(genspec))>0){
    cat("Dropping duplicated taxa: ", paste(taxalist[duplicated(genspec)], collapse=", "), "\n")
  }
  if(sum(is.na(tax$ott_id))>0){
    cat("No ott ids found for taxa: ", paste(taxalist[is.na(tax$ott_id)], collapse=", "), "\n")
  }
  tax_unique <- tax[!(duplicated(genspec) | is.na(tax$ott_id)),]
  tax_unique$ottids <- as.character(tax_unique$ott_id)
  options("scipen"=scipen, "digits"=digits)
  tax_unique[,1] <- gsub(" ", "_", tax_unique[,1])
  tax_unique[,1] <- sapply(tax_unique[,1], function(x) simpleCap(x))
  return(tax_unique)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

.nodeHeights <- function(tree, scale="YBP"){
  branchTimes <- branching.times(tree)
  TH <- max(branchTimes)
  n.tip <- length(tree$tip.label)
  branchTimes <- c(setNames(rep(0, n.tip), 1:n.tip), branchTimes)
  nH <- cbind(branchTimes[tree$edge[,1]], branchTimes[tree$edge[,2]])
  if(scale=="YBP"){
    return(nH)
  }
  if(scale=="ABS"){
    return(TH-nH)
  }
}

getBranchesSlice <- function(slice, tree){
  nH <- .nodeHeights(tree)
  present <- unname(which(apply(nH, 1, function(x) slice > x[2] & slice < x[1])))
  return(list(branches=present, nodes=tree$edge[present,2]))
}

sliceTaxonomyTable <- function(slices, tree, lookupLICAs=FALSE, ottTable=NULL){
  cache <- geiger:::.prepare.bm.univariate(tree, setNames(rep(1,length(tree$tip.label)), tree$tip.label))
  TH <- max(branching.times(tree))
  nodeSlices <- lapply(slices, function(x) getBranchesSlice(x, tree)$nodes)
  tipDescendents <- lapply(nodeSlices, function(x) cache$desc$tips[x])
  tipNames <-lapply(tipDescendents, function(x) lapply(x, function(y) cache$phy$tip.label[y]))
  if(lookupLICAs){
    tipOtts <- lapply(tipDescendents, function(x) lapply(x, function(y) ottTable$ott_id[y]))
    tipLICAs <- mclapply(tipOtts, function(x) sapply(x, function(y) rotl::taxonomy_lica(ott_ids=as.character(y))$lica$'ot:ottTaxonName'), mc.cores=2)
  }
  tipDescendents <- lapply(tipDescendents, function(x) do.call(rbind, lapply(1:length(x), function(y) cbind(y, x[[y]]))))
  if(lookupLICAs){
    tipLicaDesc <- lapply(1:length(tipDescendents), function(x) data.frame(tipDescendents[[x]],taxName=tipLICAs[[x]][tipDescendents[[x]][,1]]))
    tax <- do.call(cbind, lapply(length(slices):1, function(x) as.character(tipLicaDesc[[x]][order(tipLicaDesc[[x]][,2]),3])))
    colnames(tax) <- paste("rank",1:ncol(tax), sep="")
    tax <- data.frame(tip=cache$phy$tip.label, tax)
    } else {
    tax <- do.call(cbind, lapply(length(slices):1, function(x) tipDescendents[[x]][order(tipDescendents[[x]][,2]),1]))
    colnames(tax) <- paste("rank",ncol(tax):1, sep="")
    tax <- data.frame(tip=cache$phy$tip.label, tax)
  }
  return(tax)
}

## # Create a list of lineage tables
extractLineages <- function(otts, ncores, ottnames=NULL){
  .extractFromRaw <- function(raw){data.frame(do.call(rbind, lapply(raw[[1]]$lineage, function(x) data.frame(otName=x$name, rank=x$rank, ottid=x$ott_id, unique_name=x$unique_name) )))}
  lineages_raw <- mclapply(otts, function(x) rotl::taxonomy_taxon_info(x, include_lineage=TRUE), mc.cores=ncores)
  lineageTable <- lapply(lineages_raw, .extractFromRaw)
  if(!is.null(ottnames)){
    names(lineageTable) <- ottnames
  }
  return(lineageTable)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


matchLineages <- function(dataLineages, treeLineages, tree, ncores=1){
  matches <- parallel::mclapply(dataLineages, function(x) sapply(treeLineages, function(y) sum(match(x$ottid, y$ottid, nomatch=0)>0)), mc.cores=ncores)
  matchIds <- lapply(matches, function(x) which(x==max(x)))
  matchTaxa <- lapply(matchIds, function(x) tree$tip.label[x])
  names(matchTaxa) <- names(matchIds) <- names(matches) <- names(dataLineages)
  return(list(taxa=matchTaxa, ids=matchIds))
}

resolveDataTaxonomy <- function(matches, taxonomy){
  taxonomies <- lapply(matches$ids, function(x) taxonomy[x, ])
  .simplifyTaxonomy <- function(table, tip){
    shared <- apply(table[,-1], 2, function(x) length(unique(x))==1)
    taxrow <- table[1,-1]
    taxrow[!shared] <- NA
    cbind(tip, taxrow)
  }
  res <- lapply(1:length(taxonomies), function(x) .simplifyTaxonomy(taxonomies[[x]], names(taxonomies)[x]))
  res <- do.call(rbind, res)
  rownames(res) <- 1:nrow(res)
  res[is.na(res)] <- length((unique(unlist(taxonomy[,-1])))) + 1:length(res[is.na(res)])
  res
}

phyndrTTOL <- function(ttolObject, taxalist, timeslices, ottids=FALSE, prune=TRUE, ncores=1, infer_context=FALSE, context=NULL){
  tree     <- ttolObject$phy
  dat      <- ttolObject$dat
  lineages <- ttolObject$lineages
  tax <- sliceTaxonomyTable(timeslices, tree, lookupLICAs = FALSE)
  rm(ttolObject)
  if(!ottids){
    if(infer_context){
      ctxt <- rotl::tnrs_infer_context(names=taxalist)
      cat(paste("Inferred context: ", ctxt$context_name, sep=""))
      dataOttTable <- getOttIds(taxalist, ncores=ncores, context=ctxt$context_name)
    } else {
      dataOttTable <- getOttIds(taxalist, ncores=ncores, context=context)
    }
  } else {
    dataOttTable <- data.frame("ott_id" = taxalist)
    rownames(dataOttTable) <- 1:length(taxalist)
  }
  missing <- setdiff(1:length(taxalist), as.numeric(rownames(dataOttTable)))
  nas <- which(is.na(dataOttTable$ott_id))
  #if(ottids) otn <- NULL else otn <- taxalist[!(1:length(taxalist) %in% missing)]
  dataLineages <- extractLineages(as.character(dataOttTable$ott_id), ncores=ncores, ottnames=taxalist[!(1:length(taxalist) %in% missing)])
  matchTaxonomy <- matchLineages(dataLineages, lineages, tree, ncores=ncores)
  dataTaxonomy <- resolveDataTaxonomy(matchTaxonomy, tax)
  if(all(dataTaxonomy[1,2]==dataTaxonomy[,2])) stop("Less than 2 of the queried taxa are represented in the timetree and no dates can be obtained")
  if(prune){
    colLICA <- suppressWarnings(min(which(apply(dataTaxonomy[,-1], 2, function(x) length(unique(x))==1))))
    if(is.finite(colLICA) & colLICA != ncol(dataTaxonomy)){
      dataLICA <- unique(dataTaxonomy[, colLICA+1])
      drop <- which(!(tax[,colLICA+1] %in% dataLICA))
    ## Loop through higher taxnomic levels not in the data taxonomy and select only one species
      #if(colLICA > 2){
      #  for(j in (colLICA-1):2){
      #    drop <- c(drop, which(!(tax[,j] %in% dataTaxonomy[,j]) & duplicated(tax[,j])))
      #  }
      # }
      drop <- unique(drop)
    } else {drop <- numeric(0)}
    if(length(drop) > 0){
        ptree <- phyndr:::drop_tip(tree, tree$tip.label[drop])
        ptax <- tax[-drop, ]
    } else {
      ptree <- tree
      ptax <- tax
    }
  } else {
    ptree <- tree
    ptax <- tax
  }
  fulltax <- rbind(ptax, dataTaxonomy)
  fulltax <- fulltax[!duplicated(fulltax[,1]),]
  rownames(fulltax) <- fulltax[,1]
  fulltax <- fulltax[,-1]
  fakenames <- sapply(1:max(unique(unlist(fulltax))), function(x) paste(sample(LETTERS, 12,replace=TRUE), collapse=""))
  fulltax[1:ncol(fulltax)] <- apply(fulltax, 2, function(x) fakenames[x])
  #fulltax <- apply(fulltax,2,function(x){class(x) <- "character"; x})
  if(!is.ultrametric(ptree)){
    nH <- phytools::nodeHeights(ptree)
    exte <- which(ptree$edge[,2]<=length(ptree$tip.label))
    ptree$edge.length[exte] <- ptree$edge.length[exte] - (nH[exte,2] - min(nH[exte,2]))
  }
  phynd <- phyndr_taxonomy(ptree, taxalist[!((1:length(taxalist)) %in% nas)], fulltax)
  return(list(otts=dataOttTable, taxonomy=fulltax, phyndr=phynd, missing=taxalist[missing]))
}

getFullTree <- function(taxobj, times, ...){
  data(ttolData)
  ttPhynd <- phyndrTTOL(ttolData, taxobj, times, ...)
  newTimetree <- phyndr_sample(ttPhynd$phyndr)
  if(class(taxobj) %in% c("factor", "character")){
    synth_tree <- rotl::tol_induced_subtree(ott_ids = ttPhynd$otts$ott_id)
    contree <- congruify.phylo(newTimetree, tree, scale="PATHd8")
  }
}

makeReferenceTree <- function(tree, cores=1){
  ottTable <- getOttIds(tree$tip.label, ncores=cores)
  ottTable[,1] <- unname(sapply(gsub(" ", "_", ottTable[,1]), simpleCap))
  td <- make.treedata(tree, ottTable)
  ttolObject <- filter(td, !is.na(ott_id))
  lineages <- extractLineages(as.character(td$dat$ott_id), ncores=cores, ottnames=td$phy$tip.label)
  ttolObject$lineages <- lineages
  return(ttolObject)
}

