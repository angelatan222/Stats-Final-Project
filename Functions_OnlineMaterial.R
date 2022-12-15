## Functions for Turetsky et al. network analysis
# Please cite paper if using these functions

# Tie recriprocity
## This function calculates proportion of reciprocal ties among all ties and participants

tieReciprocated <- function(edgelist) {
  edgelist$tierecip <- NA
  for (i in 1:nrow(edgelist)){
    if(is.na(edgelist[i,]$to) | is.na(edgelist[i,]$tierecip)==FALSE) next
    match <- subset(edgelist, to==edgelist[i,]$from & from==edgelist[i,]$to)
    edgelist[i,]$tierecip <- ifelse(nrow(match)==1, 1, 0)
  }
  recipprop <- edgelist %>%
    mutate(tierecip = ifelse(to == from, NA, tierecip)) %>%
    group_by(from) %>% 
    dplyr::summarize(recip = mean(tierecip, na.rm=TRUE)) %>%
    dplyr::rename(PPID = from) %>%
    mutate_if(sapply(., is.factor), as.character)
  
  as.data.frame(recipprop)
}


# Tie recriprocity -- complete data only
## This function calculates proportion of reciprocal ties among only participants for whom we have complete data
## i.e., only among ties who have the potential to reciprocate because they were in the study 
## and completed the questionnaire

tieReciprocatedComplete <- function(edgelist) {  # function that only looks at those who completed social network survey
  completed <- unique(edgelist$from)
  edgelist$tierecip <- NA
  for (i in 1:nrow(edgelist)){
    if(is.na(edgelist[i,]$to) | is.na(edgelist[i,]$tierecip)==FALSE) next
    if(!(edgelist[i,]$to %in% completed)) next
    match <- subset(edgelist, to==edgelist[i,]$from & from==edgelist[i,]$to)
    edgelist[i,]$tierecip <- ifelse(nrow(match)==1, 1, 0)
  }
  recipprop <- edgelist %>%
    mutate(tierecip = ifelse(to == from, NA, tierecip)) %>%
    group_by(from) %>% 
    dplyr::summarize(recipcomp = mean(tierecip, na.rm=TRUE)) %>%
    dplyr::rename(PPID = from) %>%
    mutate_if(sapply(., is.factor), as.character)
  
  as.data.frame(recipprop)
}


# Function to rescale data to 0-1

rescale01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}


# Function to take social network data as collected via survey (wide format) and get clean edgelist

getEL <- function(snawide, nettype) {
  
  xlong <- snawide %>% 
    dplyr::select(c(PPID, contains(nettype))) %>%
    gather("variable","to",-PPID) %>%
    dplyr::rename(from = PPID) 
  
  el <- xlong %>%
    filter(., str_detect(variable, "Name")) %>%
    mutate(weight = filter(xlong, str_detect(variable, "Strength"))[,3]) %>%
    dplyr::select(-variable) %>%
    mutate(weight = ifelse(!is.na(to) & is.na(weight), 1, weight),
           to = ifelse(is.na(to), from, to)) %>%
    mutate(weight = ifelse(to == from, 0, weight)) %>%
    distinct()
  
  lvls <- unique(c(as.character(el$from), as.character(el$to))) 
  
  el <- el %>%
    mutate(from = factor(from, levels = lvls),
           to = factor(to, levels = lvls))
  
  el
}


# Get desired centralities in a dataframe

getCentDF <- function(graph, snawide){
  data.frame(PPID = get.vertex.attribute(graph, "name"),
             degin = degree(graph, mode = "in"),
             degout = degree(graph, mode = "out"),
             degtot = degree(graph, mode = "all"), 
             btwn = betweenness(graph, directed = T, normalized=TRUE),
             hclose = (harmonic_centrality(graph, mode = "all", weights=E(graph)$weight)/(vcount(graph)-1)),
             hclosein = (harmonic_centrality(graph, mode = "in", weights=E(graph)$weight)/(vcount(graph)-1)),
             hcloseout = (harmonic_centrality(graph, mode = "out", weights=E(graph)$weight)/(vcount(graph)-1)),
             strin = strength(graph, mode = "in", loops=FALSE),
             strout = strength(graph, mode = "out", loops=FALSE),
             strtot = strength(graph, mode = "all", loops=FALSE),
             strinstd = rescale01(strength(graph, mode = "in", loops=FALSE)),
             stroutstd = rescale01(strength(graph, mode = "out", loops=FALSE)),
             strtotstd = rescale01(strength(graph, mode = "all", loops=FALSE))) %>%
    filter(PPID %in% snawide$PPID) %>%
    mutate(strinavg = ifelse(degin==0, 0, strin/degin),
           stroutavg = ifelse(degout==0, 0, strout/degout),
           strtotavg = ifelse(degtot==0, 0, strtot/degtot)) %>%
    mutate_if(sapply(., is.factor), as.character)
}


# Grand mean center variables

gmc <- function(var){
  var-mean(var, na.rm=T)
}



