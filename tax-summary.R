#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
library(stringr)

# Recupera la informacion taxonómica
taxdata <- function(table,char=";",bind_table=F){
  cl = str_split(table[,"Consensus.Lineage"],char,n=7,T)
  colnames(cl) = c("k","p","c","o","f","g","s")
  garbage = ".__| .__"
  cl = gsub(pattern=garbage,replacement="",cl)
  if(bind_table==T){
    out = cbind(table,cl)
  } else {
    out = cl
  }
  return(out)
}

# Recupera la matriz numérica
taxmatrix <-function(table,has_otuid=T){
  if(has_otuid==T){
    out = table[,-c(1,ncol(table))]
  } else {
    out = table[,-c(ncol(table))]
  }
  return(out)
}

# Genera la sumarización de los niveles
taxsum <- function(tx,td,taxlevel){
  tab=cbind(tx,td[,taxlevel])
  colnames(tab)[ncol(tab)]<-"x"
  dfa = data.frame(matrix(ncol=ncol(tab),nrow=0),stringsAsFactors=F)
  colnames(dfa)=colnames(tab)
  for(otu in unique(tab$x)){
    sub_tab = filter(tab,x==otu)
    sum_df = colSums(sub_tab[,-ncol(sub_tab)])
    dfb = append(sum_df,otu)
    names(dfb) = colnames(tab)
    dfb = data.frame(lapply(dfb, type.convert), stringsAsFactors=FALSE)
    dfa = rbind(dfa,dfb)
    colnames(dfa)=colnames(tab)
  }
  colnames(dfa)[ncol(dfa)]<-taxlevel
  return(dfa)
}

args=commandArgs(T)
file=args[1]
table=read.table(file,header=T,sep=",",comment.char="#")

tx = taxmatrix(table)
td = taxdata(table)

for(i in c("k","p","c","o","f","g")){
  txsum = taxsum(tx,td,i)
  outfile = paste(i,file,sep=".")
  message("Generating file ",outfile)
  write.table(txsum,outfile,quote=F,sep="\t",row.names=F,col.names=T)
}
message("Done!")
