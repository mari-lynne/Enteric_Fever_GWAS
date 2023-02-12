#Volcano plot function ####

plot_vol <- function(deg, p = 0.05, FC = 0,
                     lab_type= "top", genes, top = 20,
                     title="", alpha = 0.98,
                     colours = c("#a50000","#800000","#ef5a3a","orange","yellow")) {
  
  
  #/ Checks
  if (missing(deg)){
    stop("Error: deg arg is missing. Please provide a toptable data frame")} 
  #if (!missing(genes) & !is.character(genes)){
  #stop("Error: label is not a character vector")}
  if (!missing(lab_type) & lab_type != c('top') & missing(genes)) {
    stop("Error: label requires character vector with selected genes")}
  
  ##/Set up DEG table 1)
  
  deg <- deg %>%
    mutate(reg =
             case_when(
               deg$logFC >= FC & deg$adj.P.Val <= p ~ "Sig Adj. P <0.05",
               deg$logFC <= FC & deg$adj.P.Val <= p ~ "Sig Adj. P <0.05",
               deg$logFC >= FC & deg$P.Value <= p ~ "Sig P <0.05",
               deg$logFC <= FC & deg$P.Value <= p ~ "Sig P <0.05",
               abs(deg$logFC) <= FC & deg$adj.P.Val >= p ~ "No Change",
               abs(deg$logFC) <= FC & deg$adj.P.Val <= p ~ "No Change",
               abs(deg$logFC) > FC & deg$adj.P.Val >p ~ "No Change")) %>%
    mutate(reg =
             factor(reg, levels =
                      c("Sig Adj. P <0.05", "Sig P <0.05", "No Change")))
  
  #/ 1) Define sig values
  log_p <- -log10(p)
  #Only make adj_p object if p <0.05 could do after reg table 
  
  #Adj p values #need to convert adjusted p-val to unadjusted for plot yscale
  values <- seq(0.050,0.051, by=0.00001) #vector of values approx at 0.05
  #round function not working
  deg$adj.P.Val <- round(deg$adj.P.Val, 5) #round adj.p-val col (get it out of log, so it can be filtered by values)
  adj.p <- deg[(deg$adj.P.Val %in% (values)),c("P.Value")] #Look at adj p col for matching value around 0.05
  #Maybe add an if adj.p=null set it to null so it can't be plotted, might throw an error message later tho
  adj.p <-  min(adj.p
                log_adj <- -log10(adj.p)) 
  
  #Define labels 3)
  
  if (is.null(lab_type)){#No entry for lab list, currently there's no default arg
    gene_label <- c("")
    lab_data <- NULL #Empty dataframe, need to work out how not to error
    warning("No genes highlighted")
  } else if(lab_type == "sig"){
    gene_label <- genes
    lab_data <- deg[(deg$reg == "Sig P <0.05" |deg$reg == "Sig Adj. P <0.05" )
                    & (deg$gene_name %in% genes),]
  } else if (lab_type == "ns"){
    gene_label <- genes
    lab_data <- deg[(deg$gene_name %in% genes),]}
  else if (lab_type == "top"){
    lab_data <- slice_min(deg,adj.P.Val,n=top)
    gene_label <- lab_data$gene_name
  }
  
  #Plot Volcano 4)
  vol <-
    deg %>% ggplot(aes(x=logFC,y=-log10(P.Value),label=gene_name))+
    geom_point(aes(color = P.Value, alpha=alpha))+
    labs(title = title) +
    theme_minimal() +theme(legend.position = "none") +
    geom_hline(yintercept = log_p, linetype = 2.5, alpha =0.7) +
    geom_hline(yintercept =log_adj, linetype = 2.5, alpha =0.7)+
    geom_label_repel(data=lab_data,
                     size=3.5,direction="both",nudge_y =1.6,nudge_x =0.1,
                     angle= 70,vjust= 0,segment.size= 0.5,
                     segment.color="#331002",fill="#f7f7f5")+
    scale_color_gradientn(colours = colours,
                          values=c(0,adj.p,p,1))
  return(vol)
} 


#Debugging ####  

#/ Error checking 5/9
#/ ISSUE - labels not plotting for sig settings

#replace all gene_col with gene_name
#Updates works, fix gene col setting later

#/ Error 6/9
#/ ISSUE hline sig settings not working for raw p-value on sig/ns graphs (works with top)
#Warning message:
#In min(adj.p) : no non-missing arguments to min; returning Inf

#check adj p value col in raw data might not be working
#All adj values are 0.93 - seems like an error with limma code and formation/calc in toptable

#I think the issue is, if you don't always get p-values that are less than <0.05, there is nothing to subset for the adj.p line, adj-p values are just high bc there are no sig DEGs, works with better contrasts
#Therefore if p != <0.05, plot a no adj plot, this might just have to be a separate plot chunk not sure how to do no label and a label option with one chunk. 
#Else P <0.051, plot plot with adjusted sig
#Also an issue if there just isn't close p values to 0.05, ultimately need to solve the equation
#So if we plot p values, y axis of p is done, need to get adj p=.05 into equivalent p-val (and only plot if there is adj p =<0.05 or less)




# Example:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
deg <- toptables$D0_vs_V0#toptables$h12TD_vs_h12nTD
#Modify gene names
genes_fc <- c("FCGR3A","FCGR3C","FCAR","FCGR1B","FCGR2A", "FCGR2B", "FCGR2C", "FCGR1A")
genes <- c("FCGR3A","FCGR3C","FCAR","FCGR1B","FCGR2A", "FCGR2B", "FCGR2C", "FCGR1A")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_vol(deg, lab = 'sig',genes = genes_fc)
debug(plot_vol)

