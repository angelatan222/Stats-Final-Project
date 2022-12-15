
# Data processing for Turetsky et al., Values affirmation strengthens peer social networks
# This file is run by PaperAnalyses_OnlineMaterial.Rmd

## Load libraries
install.packages("readr")
install.packages("igraph")
install.packages("reshape2")
library(tidyverse)
library(igraph)
library(CINNA)
library(reshape2)

## Load custom functions
source("Stats-Final-Project/Functions_OnlineMaterial.R")

## Read in data
att <- read.csv("Stats-Final-Project/attributes.csv") %>%  # file of participant attributes
  mutate(PPID = as.character(PPID))
sna1 <- read.csv("Stats-Final-Project/SNA_T1.csv")      # time 1 network data
sna2 <- read.csv("Stats-Final-Project/SNA_T2.csv")      # time 1 network data
prepost <- read.csv("Stats-Final-Project/PrePostMeasures.csv") %>%   # relevant survey items collected at times 1 and 2
  mutate(PPID = as.character(PPID))


## Process network data

### Friendship time 1

fel1 <- getEL(sna1, "Friend")
names(fel1)

fmat1 <- reshape2::acast(fel1, from ~ to, value.var="weight", fill=0, drop=FALSE)

fg1 <- igraph::graph.adjacency(fmat1, mode=c("directed"), weighted=TRUE, diag=FALSE)

fnetd1 <- getCentDF(fg1, sna1) %>%
  full_join(tieReciprocated(fel1), by="PPID") %>%   
  full_join(tieReciprocatedComplete(fel1), by="PPID") %>%
  mutate(Time = 1)

### Friendship time 2
  
  fel2 <- getEL(sna2, "Friend")
  
  fmat2 <- reshape2::acast(fel2, from~to, value.var="weight", fill=0, drop=FALSE)
  
  fg2 <- igraph::graph.adjacency(fmat2, mode=c("directed"), weighted=TRUE, diag=FALSE)
  
  fnetd2 <- getCentDF(fg2, sna2) %>%                
    full_join(tieReciprocated(fel2), by="PPID") %>%   
    full_join(tieReciprocatedComplete(fel2), by="PPID") %>%
    mutate(Time = 2)
  
### Study time 1

  sel1 <- getEL(sna1, "Study")
  
  smat1 <- reshape2::acast(sel1, from~to, value.var="weight", fill=0, drop=FALSE)
  
  sg1 <- igraph::graph.adjacency(smat1, mode=c("directed"), weighted=TRUE, diag=FALSE)
  
  snetd1 <- getCentDF(sg1, sna1) %>%                
    full_join(tieReciprocated(sel1), by="PPID") %>%   
    full_join(tieReciprocatedComplete(sel1), by="PPID") %>%
    mutate(Time = 1)
  
### Study time 2
  
  sel2 <- getEL(sna2, "Study")
  
  smat2 <- reshape2::acast(sel2, from~to, value.var="weight", fill=0, drop=FALSE)
  
  sg2 <- igraph::graph.adjacency(smat2, mode=c("directed"), weighted=TRUE, diag=FALSE)
  
  snetd2 <- getCentDF(sg2, sna2) %>%                
    full_join(tieReciprocated(sel2), by="PPID") %>%   
    full_join(tieReciprocatedComplete(sel2), by="PPID") %>%
    mutate(Time = 2) 
  
### Support time 1
  
  supel1 <- getEL(sna1, "Support")
  
  supmat1 <- reshape2::acast(supel1, from~to, value.var="weight", fill=0, drop=FALSE)
  
  supg1 <- igraph::graph.adjacency(supmat1, mode=c("directed"), weighted=TRUE, diag=FALSE)
  
  supnetd1 <- getCentDF(supg1, sna1) %>%                
    full_join(tieReciprocated(supel1), by="PPID") %>%   
    full_join(tieReciprocatedComplete(supel1), by="PPID") %>%
    mutate(Time = 1)
  
### Support time 2
  
  supel2 <- getEL(sna2, "Support")
  
  supmat2 <- reshape2::acast(supel2, from~to, value.var="weight", fill=0, drop=FALSE)
  
  supg2 <- igraph::graph.adjacency(supmat2, mode=c("directed"), weighted=TRUE, diag=FALSE)
  
  supnetd2 <- getCentDF(supg2, sna2) %>%                
    full_join(tieReciprocated(supel2), by="PPID") %>%   
    full_join(tieReciprocatedComplete(supel2), by="PPID") %>%
    mutate(Time = 2)

## Merge network datasets
  
### All 328 participants
  
  # wide
  all328 <- full_join(fnetd1, fnetd2, by="PPID", suffix=c(".T1.Fr",".T2.Fr")) %>%
    full_join(snetd1, by="PPID") %>%
    full_join(snetd2, by="PPID", suffix=c(".T1.St",".T2.St")) %>%
    full_join(supnetd1, by="PPID") %>%
    full_join(supnetd2, by="PPID", suffix=c(".T1.Sup",".T2.Sup")) %>%
    full_join(att, by="PPID") %>%
    full_join(prepost, by="PPID") %>%
    left_join(sna1 %>% dplyr::select(PPID, Friend1Name:Friend6Name) %>% mutate(PPID=as.character(PPID)), by="PPID") %>%
    dplyr::rename(Fr1.T1 = Friend1Name,
           Fr2.T1 = Friend2Name,
           Fr3.T1 = Friend3Name,
           Fr4.T1 = Friend4Name,
           Fr5.T1 = Friend5Name,
           Fr6.T1 = Friend6Name) %>%
    left_join(sna2 %>% dplyr::select(PPID, Friend1Name:Friend6Name) %>% mutate(PPID=as.character(PPID)), by="PPID") %>%
    dplyr::rename(Fr1.T2 = Friend1Name,
           Fr2.T2 = Friend2Name,
           Fr3.T2 = Friend3Name,
           Fr4.T2 = Friend4Name,
           Fr5.T2 = Friend5Name,
           Fr6.T2 = Friend6Name) %>%
    mutate(Recitation = as.factor(Recitation)) %>%
    mutate(RecitationC = C(Recitation, sum),        # use sum/deviation contrast coding for recitation and
           hclose.T1.Fr.gmc = gmc(hclose.T1.Fr),    # grand mean center baseline network variables so that
           btwn.T1.Fr.gmc = gmc(btwn.T1.Fr),        # intercept represents the average student
           degtot.T1.Fr.gmc = gmc(degtot.T1.Fr),    # (does not change estimates of treatment effect, just  
           degout.T1.Fr.gmc = gmc(degout.T1.Fr),    #  needed for calculation of percent change in treatment
           degin.T1.Fr.gmc = gmc(degin.T1.Fr),      #  from average student in control)
           strtotavg.T1.Fr.gmc = gmc(strtotavg.T1.Fr),
           stroutavg.T1.Fr.gmc = gmc(stroutavg.T1.Fr),
           strinavg.T1.Fr.gmc = gmc(strinavg.T1.Fr))  
    
  #long, by time
  all328long <- full_join(rbind(fnetd1, fnetd2), 
                          rbind(snetd1, snetd2), by=c("PPID","Time"), 
                          suffix=c(".Fr","")) %>%
    full_join(rbind(supnetd1, supnetd2), by=c("PPID","Time"),
              suffix=c(".St",".Sup")) %>%
    full_join(att, by="PPID") %>%
    full_join(prepost, by="PPID") %>%
    mutate(Recitation = as.factor(Recitation))

### 290 participants who completed intervention
  int290 <- subset(all328, is.na(Intervention)==FALSE)
  int290long <- subset(all328long, is.na(Intervention)==FALSE)
  
### 226 "completers" (who completed intervention and both surveys)
  complet226 <- subset(int290, CompletedPreQ==1 & CompletedPostQ==1)
  complet226long <- subset(int290long, CompletedPreQ==1 & CompletedPostQ==1)

write.csv(int290, "~/Desktop/stats 2/Final Project/int290.csv", row.names=FALSE)
write.csv(int290long, "~/Desktop/stats 2/Final Project/int290long.csv", row.names=FALSE)
write.csv(complet226, "~/Desktop/stats 2/Final Project/complet226.csv", row.names=FALSE)
write.csv(complet226long, "~/Desktop/stats 2/Final Project/complet226long.csv", row.names=FALSE)




