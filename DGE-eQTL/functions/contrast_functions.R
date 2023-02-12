#Just contrast matrix function ####
make_contrasts <- function (group, control, delim="_vs_", des_mat){
  #/ define groups and baseline to make contrasts
  
  suppressMessages(require(limma))
  
  #Checks
  if(is.null(group)) stop("Error: group arg is missing")
  
  #/ ensure unique group levels
  group <- sort(unique(as.character(group)))
  
  #Write limma code by pasting groups
  #/ if control var is present, compare all groups to control
  #/  else make all comparisons using combn function
  if (!missing(control)){
    combo <- paste0(group,"-",  control)
  } else{
    combo <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  }
  
  #/ make contrasts
  if (!missing(des_mat)){
  contrasts<- limma::makeContrasts(contrasts=combo, levels=colnames(des_mat))
  }else{
    contrasts<- limma::makeContrasts(contrasts=combo, levels=group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  }
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  #Todo: levels need to be design matrix, poss modify to have group or design option later
  #Only do lmfit step if design and fit are supplied
  
  return(contrasts)
}

#output testing
#results <- make_contrasts(group = c("V7", "V1", "D0", "V0"), control = c("V0"), des_mat = design)

#Combo with ebayes #####
contrast_2_lm <- function(group, control, delim="_vs_", design, fit, topTab="TRUE"){
  #Define groups and baseline to make contrasts
  #Input design and previous model fit model for it to work
  #Toptab gives option to output results
  
  #/ Checks
  suppressMessages(require(limma))
  if(is.null(group)) stop("Error: group arg is missing")
  
  #/ Ensure unique group levels
  group <- sort(unique(as.character(group)))
  
  #/ Define contrasts 
  if (!missing(control)){ #compare to control
    combo <- paste0(group,"-",  control)
  } else{ #make all pairwise comparisons
    combo <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])}) 
  }
  
  #/ Make contrast matrix
  if (!missing(design)){
    contrasts<- limma::makeContrasts(contrasts=combo, levels=colnames(design))
  }else{ #No design + no efit
    contrasts<- limma::makeContrasts(contrasts=combo, levels=group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  } 
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  
 
  #/ Model fit and deg results
  if (!missing(fit) & !missing(design) & topTab == "TRUE"){ #Requires fit and dmat args
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    message("Performing ebayes fit of linear model")
    top_results <- list() #TopTable Results
    for (i in colnames(contrasts)){
      top_results[[i]] <- topTable(fit2, coef = i, number = Inf)
      limma_list <- list(contrasts = contrasts, fit2 = fit2, top_results = top_results)
    }
    message("Contrast matrix (contrasts), ebayes (fit2), top DEGs (top_results) saved in list\n
          Subset list with either $ or [[]] for results")
  } else if (!missing(fit) & !missing(design) & topTab != "TRUE"){ #Without toptable option
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    message("Performing ebayes fit of linear model")
    limma_list <- list(contrasts = contrasts, fit2 = fit2)
  } else if (missing(fit) | missing(design)) {
    limma_list <- list(contrasts)
    warning("No linear model or design matrix supplied, returning contrast matrix only")
  }
  return(limma_list)
}

#Output testing 
results2 <- contrast_2_lm(group = c("V7nTD", "V1nTD", "D0nTD"), control = c("V0nTD"), fit=fit, design=design)

#6/9/22
#Error - all adj.p vals = 0.93
toptables <- results2$top_results
deg <- toptables$D0nTD_vs_V0nTD #adj.pvalues here - this has some around 0.05


results <- contrast_2_lm(group = c("h12TD"), control = c("h12nTD"),fit=fit, design = design) 
top <- results$top_results
test <- top$h12TD_vs_h12nTD


# Update 5/9/22

#/ requires option to do all comparisons automatically to control if no group is specified
#/ Add in automatic group comparison
# / not fixed yet - man input groups 4 now


make_contrasts <- function (group="all", control, delim="_vs_", des_mat){
  #/ define groups and baseline to make contrasts
  
  suppressMessages(require(limma))
  
  #Checks
  #if(is.null(group)) stop("Error: group arg is missing")
  #instead replace with deault arg 
  #Define var as all levels - control var # but if group and des is missing then it will have to error
  
  #/ ensure unique group levels
  group <- sort(unique(as.character(group)))
  
  #Write limma code by pasting groups

  if (!missing(control) & group !="all"){ #/ if control var is present, and groups specified, compare all groups to control
    combo <- paste0(group,"-",  control)
  } else if (!missing(control) & group =="all" & !missing(des_mat)){ #/ if control var is present, and no groups specified, compare all levels of design matrix as groups to control
    combo <- combn(colnames(des_mat), 2, FUN = function(x){paste0(x[1], "-", x[2])})
    combo <- colnames(design[,grep(control,colnames(des_mat))])
  } else if (missing(control)){ #make all comparisons using combn function
    combo <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])}) #This is subsetting fields to paste
  }
  
  
  #/ make contrasts
  if (!missing(des_mat)){
    contrasts<- limma::makeContrasts(contrasts=combo, levels=colnames(des_mat))
  }else{
    contrasts<- limma::makeContrasts(contrasts=combo, levels=group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  }
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  #Todo: levels need to be design matrix, poss modify to have group or design option later
  #Only do lmfit step if design and fit are supplied
  
  return(contrasts)
}

tests <- make_contrasts(control = c("time0"), des_mat = design)
