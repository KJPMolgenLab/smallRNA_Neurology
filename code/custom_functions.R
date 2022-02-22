require("kableExtra")
require("tidyverse")
require("compareGroups")
require("RColorBrewer")
require("stringr")
require("pheatmap")
require("DESeq2")
require("lm.beta")
require("viridis")


Dark8 = brewer.pal(8, "Dark2")

display_tab = function(df){
  df %>% kbl(digits = 3,align = "l") %>%
    column_spec(1, bold = T) %>%
    kable_paper(bootstrap_options = c("striped", "hover")) %>% 
    column_spec(column=2, color = spec_color(df[,1], direction=-1)) %>% 
    column_spec(column=3, color = spec_color(df[,2], direction=-1)) %>% 
    column_spec(column=4, color = spec_color(df[,3], direction=-1)) %>% 
    column_spec(column=5, color = spec_color(df[,4], direction=-1))
}


# go profiler function
getGOresults = function(geneset, genereference, organism = "mmusculus",
                        sources=c("GO:BP", "GO:MF", "GO:CC", "KEGG", "TF",
                                  "MIRNA","CORUM", "HP", "HPA")){
  require(gprofiler2)
  resgo = gost(geneset, organism =organism,
               correction_method = "gSCS",
               domain_scope = "custom",
               sources = sources,
               custom_bg = genereference,
               numeric_ns = "ENTREZGENE_ACC")
  if(length(resgo) != 0){
    return(resgo)
  } else {
    print("no significant GO terms identified")
    return(NULL)
  }
}



GOplot = function(GOtable, N, Title="GO plot"){
  if(nrow(GOtable)<N){N=nrow(GOtable)}
  GOtable = GOtable[GOtable$parents!="character(0)",]
  Tabtoplot=GOtable[order(GOtable$p_value, decreasing = F)[1:N],]
  Tabtoplot$log10pvalue=-log10(Tabtoplot$p_value)
  Tabtoplot$genperc=Tabtoplot$intersection_size/Tabtoplot$effective_domain_size

  wrapit = function(long_phrase, cutoff){
    if(nchar(long_phrase) > cutoff){
      cutpos=ceiling(str_count(long_phrase, pattern = " ")/2)
      modx = gsub(paste0("(([^ ]* ){",cutpos,"})([^ ]*)"), "\\1\n\\3", long_phrase)
      return(modx)
    } else {
      return(long_phrase)}

  }

  Tabtoplot$term_name = sapply(Tabtoplot$term_name, wrapit, cutoff=40)

  ggplot(Tabtoplot) + geom_point(aes(x =log10pvalue,
                                     y = N:1,
                                     size=precision,
                                     colour=genperc),
                                 alpha=0.7) +
    scale_colour_gradient(low="#00FF33", high ="#FF0000", guide = "colourbar")+
    labs(colour="genomic cov", size="precision")+
    xlab("- log10(p-value)") + ylab("GO term")+
    scale_size(range = c(3, 8))+
    theme_bw(base_size = 12) + ggtitle(Title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(breaks=N:1,
                       labels=Tabtoplot$term_name)
}

samesign <- function(x) {abs(sum(sign(x)))==length(x)}

multiORplot = function(datatoplot=FALSE, Pval = "Pval", Padj = "Padj", SE = "SE", beta="beta", pheno = "pheno"){
  starpval=convertpvaltostars(datatoplot[[Pval]])
  starpval[datatoplot[[Padj]]<0.05]="adj.p**"
  starpval[is.na(datatoplot[[beta]])]="n.a."
  CIUpper = datatoplot[[beta]] +1.96*datatoplot[[SE]]
  CILower = datatoplot[[beta]] -1.96*datatoplot[[SE]]
  xlim=range(c(CIUpper, CILower), na.rm=T)*1.2
  par(mar=c(5,16,5,2))
  betas = datatoplot[[beta]]

  plot(x=betas, y=1:length(betas),
       type="n", panel.first = grid(ny=NA),
       yaxt = "n", ylab="",
       xlim=xlim,
       xlab=expression(paste('log(OR)'%+-%95,"%CI")),
       main=paste(pheno))
  abline(v=0,col="black", lty=3)
  axis(2, at=1:length(betas),
       labels=base::rev(rownames(datatoplot)),
       las=1)
  arrows(x0=CILower, x1=CIUpper, y0=length(betas):1, y1=length(betas):1, col=rainbow(length(betas)), length=0, lwd=2,code = 3)
  points(y=length(betas):1, x=betas, pch=18, col="black")
  betas[is.na(betas)]=0
  text(y=(length(betas):1)+0.5, x=betas, labels=starpval, cex=0.7)
}

convertpvaltostars=function(x){
  sapply(x, function(x){ifelse(x<=0.01, "**", ifelse(x<=0.05, "*", ""))})
}

# plots eigenvalues of specified samples
EigengenePlot=function(data, Sampledata, samplesincl){
  for (i in colnames(data)){
    nf=layout(matrix(c(1:5,rep(6,5)),ncol=2),
              heights = c(12,1,1,1,1),
              widths = c(10,2))
    par(mar=c(0.2,4.1,3,1))
    barplot(data[,i], col=gsub("ME", "", i),border = NA, main=i, ylab="ME expression")
    par(mar=c(0.1,4.1,0,1))
    a=barplot(rep(1,length(data[,i])),border = NA,
              col=ann_colors[["gRNA"]][Sampledata[samplesincl,"gRNA"]], yaxt='n')
    a=barplot(rep(1,length(data[,i])),border = NA,
              col=ann_colors[["RAPA"]][Sampledata[samplesincl,"RAPA"]], yaxt='n')
    a=barplot(rep(1,length(data[,i])),border = NA,
              col=ann_colors[["DIFF"]][Sampledata[samplesincl,"DIFF"]], yaxt='n')
    a=barplot(rep(1,length(data[,i])), border = NA,
              col=ann_colors[["CellLine"]][Sampledata[samplesincl,"CellLine"]], yaxt='n')
    par(mar=c(0,0,0,0))
    plot.new()
    legend(0,0.5, legend = names(ann_colors[["gRNA"]]),
           fill = ann_colors[["gRNA"]], xpd=T,
           bty = "n")
    legend(0,0.375, legend = names(ann_colors[["RAPA"]]),fill = ann_colors[["RAPA"]], xpd=T,bty = "n")
    legend(0,0.25, legend = names(ann_colors[["DIFF"]]),fill = ann_colors[["DIFF"]], xpd=T,bty = "n")
    legend(0,0.125, legend = names(ann_colors[["CellLine"]]),fill = ann_colors[["CellLine"]], xpd=T,bty = "n")
  }
}


lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

