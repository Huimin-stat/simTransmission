#' @import stats 
#' @import phybase
#' @import network
#' @import GGally
#' @import RColorBrewer
#' @import graphics
#' @import sna
#' @importFrom dplyr if_else
NULL

#' Generate random numbers from uniform distribution with input parameters in vector format
#'
#' @param input a vector containing 3 elements: amount of random numbers to be generated, lower bound and upper bound of uniform distribution
#' @return a vector of generated random numbers
Runif <- function(input){
  c <- runif(input[1], min = input[2], max = input[3])
  c
}

#' Generate transmission information of next generation with consideration of asyptomatic/unsampled individuals
#' 
#' @param currentmatrix information matrix of current generation 
#' @param inf_rate per day infective rate per infected capita; could be a constant or vector
#' @param inf_rate_timeInt time interval correspond to the start and end time of each infectious rate; matrix
#' @param diag_rate per day diagnosis rate per capita for symptomatic/sampled individuals; constant
#' @param rec_rate_asy per day recovery rate per capita for asymptomatic/unsampled individulas; constant
#' @param asy_prop aysmptomatic propotion among infected individuals; constant
#' @param rec_rate_sym per day recovery rate per capita for symptomatic/unsampled individulas after diagnosis; constant
#' @param death_prop proportion of dead cases among diagnosed people; constant
#' @param death_rate death rate per day per capita (inverse of the average time/days to death); constant
#' @param inf_rate_diag infection rate per day per capita for symptomatic individuals after diagnosis; constant
#' @param unif_input a n*3 matrix containing input parameters to generate random numbers/times from uniform distribution
#' @param ntrans number of infected individuals from the begining of outbreak (first generation) to the current generation
#' @return a list containing transmission information matrix and input parameter matrix to generate random numbers/times of next generation from uniform distribution
transNext = function(currentmatrix, inf_rate, inf_rate_timeInt, diag_rate, rec_rate_asy, rec_rate_sym, death_prop, death_rate, inf_rate_diag, asy_prop, unif_input, ntrans)
{
  nrow = dim(currentmatrix)[1]
  m = nrow(inf_rate_timeInt)
  
  index = 0
  nextmatrix = matrix(0,1,8)
  nextunif_input = matrix(0,1,3)
  
  for(i in 1:nrow){
    if(currentmatrix[i,3] == 0) next
    inft = sort(unlist(apply(unif_input[(1+(i-1)*(m+1)):(i*(m+1)),],1,Runif)))
    asy_prop_v = rep(asy_prop, currentmatrix[i,3])
    death_prop_v = rep(death_prop, currentmatrix[i,3])
    asy = if_else(asy_prop_v <= runif(currentmatrix[i,3]), 0, 1)
    
    dead = if_else(asy == rep(0,currentmatrix[i,3]),if_else(death_prop_v <= runif(currentmatrix[i,3]), 0, 1),0)
    diagt = if_else(asy == rep(0,currentmatrix[i,3]), rexp(currentmatrix[i,3], rate=diag_rate), NaN)
    rect = rexp(currentmatrix[i,3], rate=(1-dead)*(rec_rate_sym*(1-asy)+rec_rate_asy*asy)+dead*death_rate)
    remt = if_else(asy == rep(0,currentmatrix[i,3]), rect+diagt, rect)
    m = nrow(inf_rate_timeInt)
    rect_v = rep(if_else(asy == rep(0,currentmatrix[i,3]), diagt+inft, rect+inft), each=m)
    inft_v = rep(inft, each=m)
    inf_rate_timeInt_v = matrix(rep(t(inf_rate_timeInt),currentmatrix[i,3]), ncol=ncol(inf_rate_timeInt), byrow=TRUE)
    inf_index = if_else((inft_v >= inf_rate_timeInt_v[,2] | rect_v <= inf_rate_timeInt_v[,1]), 0, 1)
    
    ninf_v = rbind(matrix(rpois(m*currentmatrix[i,3],
                                lambda = inf_rate*inf_index*(pmin(rect_v,inf_rate_timeInt_v[,2])-pmax(inft_v,inf_rate_timeInt_v[,1]))), nrow = m), 
                   rpois(currentmatrix[i,3],lambda = (1-asy)*rep(inf_rate_diag,currentmatrix[i,3])*rect))
    n = colSums(ninf_v)
    
    unif_lp = as.vector(rbind(matrix(inf_index*pmax(inft_v,inf_rate_timeInt_v[,1]), nrow = m), if_else(asy == 0,(diagt+inft),0)))
    unif_up = as.vector(rbind(matrix(inf_index*pmin(rect_v,inf_rate_timeInt_v[,2]), nrow = m), (1-asy)*(remt+inft)))
    temunif_input = matrix(c(as.vector(ninf_v),unif_lp, unif_up), ncol = 3)
    
    x = cbind(inft,inft+remt, n, ntrans+1:currentmatrix[i,3], currentmatrix[i,4], asy, inft+diagt, dead)
    nextmatrix = rbind(nextmatrix,x)
    nextunif_input = rbind(nextunif_input,temunif_input)
    ntrans = ntrans + currentmatrix[i,3]
  }
  
  nextmatrix = nextmatrix[-1,]
  nextunif_input = nextunif_input[-1,]
  
  results = list(nextmatrix,nextunif_input)
  names(results) = c("nextmatrix", "nextunif_input")
  return(results)
}

#' Generate the first infected individual in a transmission tree with consideration of asyptomatic/unsampled individuals
#' 
#' @param inftime the start time of infection/outbreak
#' @param inf_rate per day infective rate per infected capita; could be a constant or vector
#' @param inf_rate_timeInt time intervals correspond to the start and end time of each infectious rate; matrix
#' @param diag_rate per day diagnosis rate per capita for symptomatic/sampled individuals; constant
#' @param rec_rate_asy per day recovery rate per capita for asymptomatic/unsampled individulas; constant
#' @param rec_rate_sym per day recovery rate per capita for symptomatic/unsampled individulas after diagnosis; constant
#' @param death_prop proportion of dead cases among diagnosed people; constant
#' @param death_rate death rate per day per capita (inverse of the average time/days to death); constant
#' @param inf_rate_diag infection rate per day per capita for symptomatic individuals after diagnosis; constant
#' @param asy_prop aysmptomatic propotion among infected individuals; constant
#' @param currentnode current node number of infected individual
#' @param ancestor node number of individual who infected the current person
#' @return a list containing transmission information vector and input parameter vector to generate random numbers/times of next generation from uniform distribution
transCurrent = function(inftime=0, inf_rate, inf_rate_timeInt, diag_rate, rec_rate_asy, rec_rate_sym, death_prop, death_rate, inf_rate_diag, asy_prop, currentnode, ancestor)
{
  # asy = 0 indicating the infected individual is sympoatic and asy = 1 implying asympomatic
  asy = if_else(asy_prop <= runif(1), 0, 1)
  dead = if_else(asy == 0,if_else(death_prop <= runif(1), 0, 1),0)
  diagt = if_else(asy == 0, rexp(1, rate=diag_rate), NaN)
  rect = rexp(1,rate=(1-dead)*(rec_rate_sym*(1-asy)+rec_rate_asy*asy)+dead*death_rate)
  remt = if_else(asy == 0, rect+diagt, rect)
  m = nrow(inf_rate_timeInt)
  rect_v = rep(if_else(asy == 0, diagt+inftime, rect+inftime), m)
  inft_v = rep(inftime, m)
  inf_index = if_else((inft_v >= inf_rate_timeInt[,2] | rect_v <= inf_rate_timeInt[,1]), 0, 1)
  
  ninf_v = c(rpois(m,lambda = inf_rate*inf_index*(pmin(rect_v,inf_rate_timeInt[,2])-pmax(inft_v,inf_rate_timeInt[,1]))),rpois(1,lambda = (1-asy)*inf_rate_diag*rect))
  ninf = sum(ninf_v)
  unif_lp = c(inf_index*pmax(inft_v,inf_rate_timeInt[,1]), if_else(asy == 0,(diagt+inftime),0))
  unif_up = c(inf_index*pmin(rect_v,inf_rate_timeInt[,2]), (1-asy)*(remt+inftime))
  unif_input = matrix(c(ninf_v,unif_lp,unif_up), ncol = 3)
  
  transnet = c(inftime,inftime+remt,ninf, currentnode, ancestor, asy, inftime+diagt, dead)
  results = list(transnet,unif_input)
  names(results) = c("transnet", "unif_input")
  return(results)
}


#' Simulate a transmission tree considering asymptomatic/unsampled individulas 
#' 
#' Takes in infection rate, recovery rate of symptomatic and asymptomatic group, asymptomatic proportion among infeted people and number of generation to generate a transmission tree
#' 
#' @param inf_rate per day infective rate per infected capita
#' @param inf_rate per day infective rate per infected capita; could be a constant or vector
#' @param inf_rate_time time points corresponds to the start time of each infectious rate when it is not a constant
#' @param diag_rate per day diagnosis rate per capita for symptomatic/sampled individuals; constant
#' @param rec_rate_asy per day recovery rate per capita for asymptomatic/unsampled individulas; constant
#' @param rec_rate_sym per day recovery rate per capita for symptomatic/unsampled individulas after diagnosis; constant
#' @param death_prop proportion of dead cases among diagnosed people; constant
#' @param death_rate death rate per day per capita (inverse of the average time/days to death); constant
#' @param inf_rate_diag infection rate per day per capita for symptomatic individuals after diagnosis; constant
#' @param asy_prop aysmptomatic propotion among infected individuals; constant
#' @param ngeneration number of generation that get infected counted from the first infected individual
#' @return transmission matrix including infection time, recovery/dead time, number of people infected by current individual, current node number of current individual, node number of individual who infected the current person (ancestor), asymptomatic information, diagnosed time and death information of current person
#' @export 
simTransmission = function(inf_rate, inf_rate_time = 0, diag_rate, rec_rate_asy, rec_rate_sym, death_prop, death_rate, inf_rate_diag, asy_prop, ngeneration)
{
  if (length(inf_rate) != length(inf_rate_time)) {
    stop("infection rates and their start time points have different length; check the input 'inf_rate' and 'inf_rate_time'.")
  }
  
  # convert the infetion rate changing time into intervals
  m = length(inf_rate_time)
  inf_rate_time[m+1] = Inf
  inf_rate_timeInt = matrix(0,m,2)
  for (i in 1:m) {
    inf_rate_timeInt[i,] = c(inf_rate_time[i], inf_rate_time[i+1])
  }
  
  #initialize transmatrix with the first infected individual of infection time 0
  transmatrix = matrix(0,1,8)
  colnames(transmatrix) = c("inftime","rectime/deadtime","ninf","current","ancestor","asymptomatic","diagtime","dead")
  currentout = transCurrent(inftime=0, inf_rate, inf_rate_timeInt, diag_rate, rec_rate_asy, rec_rate_sym, death_prop, death_rate, inf_rate_diag,asy_prop, currentnode=1, ancestor=-1)
  unif_input = currentout$unif_input
  transmatrix[1,] = currentout$transnet
  
  currentmatrix = transmatrix
  ngen = 1
  
  for(i in 1:ngeneration){
    if(length(dim(currentmatrix)) == 0) currentmatrix = t(as.matrix(currentmatrix))
    if(sum(currentmatrix[,3]) == 0) break;
    nextout = transNext(currentmatrix, inf_rate, inf_rate_timeInt, diag_rate, rec_rate_asy, rec_rate_sym, death_prop, death_rate, inf_rate_diag, asy_prop, unif_input, ntrans = dim(transmatrix)[1])
    nextm = nextout$nextmatrix
    unif_input = nextout$nextunif_input
    if(i == ngeneration) nextm[,3] = 0
    transmatrix = rbind(transmatrix,nextm)
    currentmatrix = nextm
    ngen = ngen + 1
  }
  
  row.names(transmatrix) = rep("",dim(transmatrix)[1])
  
  if(nrow(transmatrix) == 1) {warning('Inital infective individual recovered without infecting anyone. Simulate another transmission tree.')}
  else if (ngen < ngeneration) {warning(paste0('Pandemic is over within ', paste0(ngeneration, paste0(' generation, it only survive for ', paste0(ngen, ' generations. Simulate another transmission tree')))))}
  else {return(transmatrix)}
}


#' Generate summary table and barplots (optional) of asymptomatic & symptomatic cases at given times since pandemic starts
#' 
#' @param transmatrix information matrix of transmission including infection time, recovery/dead time, number of people infected by current individual, current node number of current individual, node number of individual who infected the current person (ancestor), asymptomatic information, diagnosed time and death information of current person
#' @param times a vector of truncated times (in days) during pandemic 
#' @param plot a logical parameter indicating whether to produce barplot or not (It's recommended to set to FALSE if vector "times" has a length longer than 15)
#' @return barplots (optional) and summary tables illustrating number of cases of six groups (asymptomatic_infectious, asymptomatic_recovered, symptomatic_latent, symptomatic_diagnosed, symptomatic_recovered, dead) by truncated times.
#' @export
cases.plot = function(transmatrix, times, plot = FALSE)
{
  if (max(transmatrix[,2]) <= min(times)) 
  {
    cases <- c(0, sum(transmatrix[,6]), 0, nrow(transmatrix)-sum(transmatrix[,6]), nrow(transmatrix)-sum(transmatrix[,6])-sum(transmatrix[,8], transmatrix[,8]))
    xx <- barplot(cases, space = c(10,0,0,0), col = brewer.pal(n = 4, name = "Set2"), beside = TRUE,
                  names.arg = c("Asy_Infectious", "Asy_Recovered", "Sym_Latent", "Sym_Diagnosed", "Sym_Recovered", "Dead"),
                  main = paste0("Asymptomatic & Symptomatic cases at and after ", paste0(min(times), " days")),
                  xlab = "Groups of infected individuals", ylab = "number of cases", 
                  ylim = c(0, 1.1*max(cases)))
    text(xx, cases, label = cases, cex = 1, pos = 3, col = "red")
    stop('Pandemic is over before any of the given time: all the symptomatic and asymptomatic cases have recovered or been dead')
  }
  
  n = length(times)  
  times = times[order(times)]
  n_asy_infectious = numeric(n)
  n_asy_recoverd = numeric(n)
  n_sym_latent = numeric(n)
  n_sym_confirmed = numeric(n)
  n_sym_recovered = numeric(n)
  n_dead = numeric(n)
  for (i in 1:n)
  {
    index = numeric()
    index1 = numeric()
    index = which(transmatrix[,1]>times[i])
    if (length(index)==0) {result = transmatrix
    } else {result = transmatrix[-index,]}
    if (is.matrix(result) == FALSE) result = t(as.matrix(result))
    result1 = result[which(result[,2]>times[i]),]
    result2 = result[which(result[,2]<=times[i]),]
    if (is.matrix(result1) == FALSE) result1 = t(as.matrix(result1))
    if (is.matrix(result2) == FALSE) result2 = t(as.matrix(result2))
    n_asy_infectious[i] = length(which(result1[,6]>0))
    n_asy_recoverd[i] = length(which(result[,6]>0))-n_asy_infectious[i]
    n_sym_latent[i] = length(which(((result1[,6] == 0)) & (result1[,7] > times[i])))
    n_sym_confirmed[i] = length(which(result[,6] == 0))-n_sym_latent[i] 
    n_sym_recovered[i] = length(which((result2[,6] == 0) & (result2[,8] == 0)))
    n_dead[i] = length(which(result2[,6] == 0))-n_sym_recovered[i]
  }
  cases = rbind(n_asy_infectious,n_asy_recoverd,n_sym_latent,n_sym_confirmed,n_sym_recovered,n_dead)
  colnames(cases) = c(paste0(times, " days"))
  
  if (plot==TRUE) {
    barplot(cases, xlab = "Time since outbreak first started", col = brewer.pal(n = 6, name = "Set2"),
            main = "Asymptomatic & Symptomatic cases at given time since outbreak",
            ylab = "number of cases", beside = TRUE, cex.names = 0.8)
    legend("topleft", c("Asy_infectious", "Asy_recovered", "Sym_latent", "Sym_diagnosed", "Sym_recovered", "Dead"), fill = brewer.pal(n = 6, name = "Set2"), cex = 0.8)
  }
  return(cases)
}


#' Truncate transmission matrix by a given time and produce total number and index of cases belong to each category.
#' 
#' Takes in transmission matrix and trucated it at given time. Then generate the truncated transmission matrix and latent/infectious and recovered cases for both syptomatic and asymptomatic groups, as well as diagnosed and dead cases by a truncated/given time
#' @param transmatrix information matrix of transmission including infection time, recovery/dead time, number of people infected by current individual, current node number of current individual, node number of individual who infected the current person (ancestor), asymptomatic information, diagnosed time and death information of current person
#' @param time a truncated time  
#' @return a time truncated transmission matrix along with total number and index of Symptomatic_latent, Symptomatic_diagnosed, Symptomatic_recovered, Asymptomatic_infectious, Asymptomatic_recovered and dead cases at the truncated time
#' @export
timeTruncate = function(transmatrix, time)
{
  if (max(transmatrix[,2]) <= time) 
  {
    print(paste0("The pandemic is over within the truncated time: ", paste0(time, " days.")))
    index.sym = which((transmatrix[,6]==0) & (transmatrix[,8]==0))
    index.asy = which(transmatrix[,6]>0)
    index.dead = which(transmatrix[,8]>0)
    n_sym = length(index.sym)
    n_asy = length(index.asy)
    n_dead = length(index.dead)
    results = list(transmatrix,n_sym,n_sym,n_asy,n_dead,index.sym, index.sym, index.asy, index.dead)
    names(results) = c("truncTM","n_sym_diagnosed", "n_sym_recovered","n_asy_recovered","n_dead","index_sym_diagnosed","index_sym_recovered", "index_asy_recovered","index_dead")
    return(results)
  }
  
  index = which(transmatrix[,1]>time)
  if (length(index)==0) {result = transmatrix
  } else {result = transmatrix[-index,]}
  
  
  if (length(index)!=0 && min(index)<=dim(result)[1]){
    for(i in min(index):dim(result)[1]){
      result[i,5] = which(result[,4]==result[i,5])
    }
    result[,4] = 1:dim(result)[1]
    x = table(result[-1,5])
    index = as.numeric(names(x))
    result[index,3] = x
    result[-index,3] = 0
  }
  
  index.sym.unc = which(result[,6]==0 & result[,7]>time)
  index.sym.con = which(result[,6]==0 & result[,7]<=time)
  index.sym.rec = which(result[,6]==0 & result[,2]<=time & result[,8]==0)
  index.asy.unc = which(result[,6]>0 & result[,2]>time)
  index.asy.con = which(result[,6]>0 & result[,2]<=time)
  index.dead = which(result[,2]<=time & result[,8]>0)
  result[which(result[,2]>time),2] = time
  result[which(result[,7]>time),7] = time
  n_sym_unconfirm = length(index.sym.unc)
  n_sym_confirmed = length(index.sym.con)
  n_sym_recovered = length(index.sym.rec)
  n_asy_infectious = length(index.asy.unc)
  n_asy_recovered = length(index.asy.con)
  n_dead = length(index.dead)
  
  results = list(result, n_sym_unconfirm, n_sym_confirmed, n_sym_recovered, n_asy_infectious, n_asy_recovered, n_dead,index.sym.unc, index.sym.con, index.sym.rec, index.asy.unc, index.asy.con, index.dead)
  names(results) = c("truncTM", "n_sym_latent", "n_sym_diagnosed", "n_sym_recovered","n_asy_infectious", "n_asy_recovered", "index_dead","index_sym_latent", "index_sym_diagnosed","index_sym_recovered", "index_asy_infectious", "index_asy_recovered", "index_dead")
  return(results)
}


#' Generate transmission matrix that can be used to plot transmission network
#' 
#' Takes in transmission matrix and convert it to a matrix giving network structure
#' 
#' @param transmatrix information matrix of transmission including infection time, recovery/dead time, number of people infected by current individual, current node number of current individual, node number of individual who infected the current person (ancestor), asymptomatic information, diagnosed time and death information of current person
#' @return network structure matrix
#' @export 
toAdj = function(transmatrix)
{
  adjmatrix = matrix(0, dim(transmatrix)[1], dim(transmatrix)[1])
  for(i in 2:dim(transmatrix)[1]) adjmatrix[transmatrix[i,5],i] = 1
  adjmatrix
}


#' Plot transmission network
#' 
#' Takes in transmission matrix and plot a transmission network
#' @param transmatrix information matrix of transmission including infection time, recovery/dead time, number of people infected by current individual, current node number of current individual, node number of individual who infected the current person (ancestor), asymptomatic information, diagnosed time and death information of current person
#' @param index_sym_latent index of symptomatic individuals who are in latent period
#' @param index_sym_diagnosed index of symptomatic individuals who have been diagnosed
#' @param index_sym_recovered index of symptomatic individuals who have recovered 
#' @param index_asy_infectious index of asymptomatic individuals who are still infectious (not recovered).
#' @param index_asy_recovered index of asymptomatic individuals who have recovered
#' @param index_dead index of individuals who was dead from the pandemic
#' @return a plot of transmission network
#' @export
plotTransNet=function(transmatrix, index_sym_latent = integer(0), index_sym_diagnosed, index_sym_recovered, index_asy_infectious = integer(0), index_asy_recovered, index_dead)
{
  adjmatrix = toAdj(transmatrix)
  tran_network = network(adjmatrix)
  col = rep("Sym_Diagnosed_Infectious", dim(transmatrix)[1])
  col[index_asy_infectious] = "Asy_Infectious"
  col[index_asy_recovered] = "Asy_Recovered"
  col[index_sym_latent] = "Sym_Latent"
  col[index_sym_recovered] = "Sym_Diagnosed_Recovered"
  col[index_dead] = "Dead"
  
  ggnet2(tran_network, color = col, size.legend ="ninf", shape = c(17 ,rep(19 ,dim(transmatrix)[1]-1)),
         size=transmatrix[,3], node.label = transmatrix[,4], color.palette = "Set2",
         label.size = 3, legend.position = "bottom", alpha = 0.8, edge.alpha = 0.5)
}


#' Find the infected individual and its next generation
#' 
#' Takes in transmission matrix, node of infected individual and return the information matrix of input node of infected individual as well as the individuals directly infected by it.
#' @param node node number of infected individual
#' @param transmatrix information matrix of transmission 
#' @return a transmission matrix including infection time, removal time and node number of current individual and the individuals infected by the input node.
findInfections = function(node, transmatrix)
{
  n = transmatrix[node,3] + 1
  infections = matrix(0,n,3)
  colnames(infections) = c("inftime","rectime","nodes")
  infections[1,3] = node
  infections[2:n,3] = which(transmatrix[,5] == node)
  infections[,1] = transmatrix[infections[,3],1]
  infections[,2] = transmatrix[infections[,3],2]
  
  infections
}

#' Simulate sub-transmission tree from one ancestor to its next generation
#' @param inf_time a vector of infection times
#' @param rec_time a vector of recovery times
#' @param nodes nodes number of infected individuals
#' @return sub-transmission tree in phylip format
subTree = function(inf_time, rec_time, nodes)
{
  n = length(nodes)
  
  if(n == 1) paste(nodes[1],":",rec_time[1]-inf_time[1],sep="")
  else if(n == 2){
    x = paste("(",nodes[1],":",rec_time[1]-inf_time[n],",",nodes[n],":",rec_time[n]-inf_time[n],")",sep="")
    x = paste(x,":",inf_time[2]-inf_time[1],sep="")
  } 
  else{
    x = paste("(",nodes[1],":",rec_time[1]-inf_time[n],",",nodes[n],":",rec_time[n]-inf_time[n],")",sep="")
    for(i in (n-1):2){
      x = paste("(",x,":", inf_time[i+1]-inf_time[i],",",nodes[i],":",rec_time[i]-inf_time[i],")",sep="")
    }
    x = paste(x,":",inf_time[2]-inf_time[1],sep="")
  }
  
  x
}

#' Generate transmission tree in phylip format
#' 
#' @param transmatrix information matrix of transmission 
#' @return transmission tree in phylip format
#' @export
modelTreefromTransmission = function(transmatrix)
{
  nnodes = dim(transmatrix)[1]
  x = findInfections(1,transmatrix)
  tree = subTree(x[,1],x[,2],x[,3])
  
  for(i in 2:nnodes){
    if(transmatrix[i,3] == 0) next
    x = findInfections(i,transmatrix)
    tree1 = subTree(x[,1],x[,2],x[,3])
    tree = gsub(paste(i,"\\:[0-9]+.[0-9]+",sep=""),tree1,tree)
  }
  tree = paste(tree,";",sep="")
  tree
}


#' Calculate the distance/time of input node from the first infected individual
#' 
#' @param node a number of node
#' @param nodematrix a matrix that describes the relationships among nodes and corresponding branch lengths
#' @return distance/time of input node from the first infected individual
distNodeToroot = function(node, nodematrix)
{
  distance = 0
  father = node
  
  while(nodematrix[father,1] > 0){
    distance = distance + nodematrix[father,4]
    father = nodematrix[father,1]
  }
  distance
}


#' Simulate phylogenetic/coalescent tree from transmission tree
#' 
#' @param ntree number of coalescent tree that one wants to generate
#' @param modeltree transmission tree in phylip format
#' @param theta product of effective sample size and mutation rate
#' @return coalescent tree in phylip format
#' @export
simCoalModeltree = function (ntree, modeltree, theta) 
{
  nspecies = length(species.name(modeltree))
  x = read.tree.nodes(modeltree, name=1:nspecies)
  nodematrix = x$nodes
  name = x$names
  rootnode = which(nodematrix[,1] == -9)
  seq = rep(1,nspecies)
  nodematrix[,5] = theta
  
  distance = 1:nspecies
  for(i in 1:nspecies){
    distance[i] = distNodeToroot(i,nodematrix)
  }
  
  maxdist = max(distance)
  
  #make a clocktree
  nodematrix[1:nspecies,4] = nodematrix[1:nspecies,4] + maxdist - distance
  
  #simulate coalescent trees
  coaltree = 1:ntree
  for(i in 1:ntree){
    x = sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,name)
    
    newtree = read.tree.nodes(x$gt, name = 1:nspecies)$nodes
    newtree[1:nspecies,4] = newtree[1:nspecies,4] - (maxdist - distance)
    
    coaltree[i] = write.subtree(rootnode, newtree, taxaname = name, root = rootnode)
  }
  
  coaltree
}


#' Plot phylogenetic/transmission tree
#' 
#' @param tree tree in phylip format
#' @return phylogenetic tree plot
#' @export
plotTree <- function(tree)
{
  tree.plot(tree)
}

