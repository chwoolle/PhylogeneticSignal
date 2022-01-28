#Measuring Phylogenetic Signal in Fossil Data

#Authors: Hank Woolley, Jeff Thompson, Becky Wu

#Citation: Woolley et al. (In Review) A biased fossil record can preserve reliable phylogenetic signal

#=======================================================================================
#CONTENTS

#1. Consistency Index in morphological characters in fossil data

#2. Retention Index in morphological characters in fossil data

#3. Delta-Statistic in morphological characters in fossil data

#=======================================================================================

#1. Consistency Index in morphological characters in fossil data

# Below, we have written code for calculating the Consistency Index of individual characters in a given morphological character dataset.

# Load required libraries
library(ape)
library(phangorn)
library(TreeSearch)
library(phytools)

#Set working directory to where data files are held 
setwd("~/Directory")

#Import trees from directory (.tre file) into R. NOTE: Fossil taxa need to be removed from the .tre file and the character matrix. We used Mequite to do so, but this can also be achieved using the droptip function in R.
Trees_Drop=read.nexus(".tre file") 

#The below line converts the tree to a multiphylo object
class(Trees_Drop)="multiPhylo"

#Import the character matrix (.nex file), with the fossil taxa dropped previously in Mesquite, as a .phydat object

Matrix_Phydat<-ReadAsPhyDat(".nex file")

#For loop to put CI of n characters of i Most Parsimonious Trees into a table/data.frame
Result1=NULL
for (i in 1:(number of MPTs)) {
  Result1[[i]] = CI(Trees_Drop[[i]], Matrix_Phydat, sitewise=TRUE)
  
}
Result1=as.data.frame(Result1)

#Transpose the table so columns=Trees
Result1=t(Result1)

#Name Columns and Rows
N=length(Result1[1,])
colnames(Result1) = c(paste("Tree_", 1:N, sep = ""))
row.names(Result1) = c(paste("Char_", 1:length(Result1[,1]), sep=""))

#Save as csv file.
write.csv(Result1, "CI_Results.csv")

#-------------------------------------------------------------------------------------------
#Below, we separate the characters into bins to test for differences between 
#more frequently preserved fossil data and less frequently preserved fossil data

CI=read.csv("CI_Results.csv", row.names = 1)

#Example partition of morphological characters:

Premaxilla=CI[1:16, ]
Cranial=CI[c(17:103, 118:195, 251:331), ]
Maxilla=CI[104:117, ]
Vomer=CI[196:213, ]
Palatine=CI[214:239, ]
Pterygoid=CI[240:250, ]
Mandibular=CI[332:387, ]
Dental=CI[388:411, ]
Tooth_bearing=c(Premaxilla, Maxilla, Vomer, Palatine, Pterygoid, Mandibular, Dental)
Dermal=CI[534:548,]
Other_Ossifications=CI[c(412:427, 450, 549:572), ]

Overrep=c(Tooth_bearing, Axial)
Underrep=c(Cranial, Pectoral, Pelvic, Appendicular, Dermal, Other_Ossifications)

#Create boxplot of CI values for anatomical partitions of data
boxplot(Overrep, Underrep, ylim=c(0, 1.0))

#-------------------------------------------------------------------------------------------------
#NON-PARAMETRIC STATISTICAL TESTS

#Mann Whitney aka Wilcoxon Rank-Sum Test
#Ho: The median CI values of Overrepresented anatomical elements in the fossil record are not statistically different from CI values of Underrepresented anatomical elements in the fossil record.

#Ha: The median CI values of Overrepresented anatomical elements in the fossil record are statistically different from the median CI values of Underrepresented anatomical elements in the fossil record.

#two-sided test
?wilcox.test
wilcox.test(Overrep, Underrep, mu=0, alt = "two.sided", conf.int = TRUE, conf.level = 0.95, paired = FALSE, exact = FALSE, correct = TRUE)

#Kolmogorov-Smirnov Test:

#Ho: The CI values of Overrepresented anatomical elements in the fossil record and the CI values of Underrepresented anatomical elements in the fossil record are drawn from populations of data with the same cumulative distribution.

#Ha: The CI values of overrepresented anatomical elements in the fossil record and the CI values of underreprensented anatomical elements in the fossil record are drawn from populations of data that have different cumulative distributions. 


ks.test(Overrep, Underrep)

#=======================================================================================

#2. Retention Index in morphological characters in the fossil record

# Below, we have written code for calculating the Retention Index of individual characters in a given morphological character dataset.

# Load required libraries
library(ape)
library(phangorn)
library(TreeSearch)
library(phytools)

#Set working directory to where data files are held 
setwd("~/Directory")

#Import trees from directory (.tre file) into R. NOTE: Fossil taxa need to be removed from the .tre file and the character matrix. We used Mequite to do so, but this can also be achieved using the droptip function in R.
Trees_Drop=read.nexus(".tre file") 

#The below line converts the tree to a multiphylo object
class(Trees_Drop)="multiPhylo"

#Import the character matrix (.nex file), with the fossil taxa dropped previously in Mesquite, as a .phydat object

Matrix_Phydat<-ReadAsPhyDat(".nex file")

#For loop to put RI of n characters of i Most Parsimonious Trees into a table/data.frame
Result1=NULL
for (i in 1:(number of MPTs)) {
  Result1[[i]] = RI(Trees_Drop[[i]], Matrix_Phydat, sitewise=TRUE)
  
}
Result1=as.data.frame(Result1)

#Transpose the table so columns=Trees
Result1=t(Result1)

#Name Columns and Rows
N=length(Result1[1,])
colnames(Result1) = c(paste("Tree_", 1:N, sep = ""))
row.names(Result1) = c(paste("Char_", 1:length(Result1[,1]), sep=""))

#Save as csv file.
write.csv(Result1, "RI_Results.csv")

#-------------------------------------------------------------------------------------------
#Below, we separate the characters into bins to test for differences between 
#more frequently preserved fossil data and less frequently preserved fossil data

RI=read.csv("RI_Results.csv", row.names = 1)

#Example partition of morphological characters:

Premaxilla=RI[1:16, ]
Cranial=RI[c(17:103, 118:195, 251:331), ]
Maxilla=RI[104:117, ]
Vomer=RI[196:213, ]
Palatine=RI[214:239, ]
Pterygoid=RI[240:250, ]
Mandibular=RI[332:387, ]
Dental=RI[388:411, ]
Tooth_bearing=c(Premaxilla, Maxilla, Vomer, Palatine, Pterygoid, Mandibular, Dental)
Dermal=RI[534:548,]
Other_Ossifications=RI[c(412:427, 450, 549:572), ]

Overrep=c(Tooth_bearing, Axial)
Underrep=c(Cranial, Pectoral, Pelvic, Appendicular, Dermal, Other_Ossifications)

#Create boxplot of RI values for anatomical partitions of data
boxplot(Overrep, Underrep, ylim=c(0, 1.0))

#-------------------------------------------------------------------------------------------------
#NON-PARAMETRIC STATISTICAL TESTS

#Mann Whitney aka Wilcoxon Rank-Sum Test
#Ho: The median RI values of Overrepresented anatomical elements in the fossil record are not statistically different from RI values of Underrepresented anatomical elements in the fossil record.

#Ha: The median RI values of Overrepresented anatomical elements in the fossil record are statistically different from the median RI values of Underrepresented anatomical elements in the fossil record.

#two-sided test
?wilcox.test
wilcox.test(Overrep, Underrep, mu=0, alt = "two.sided", conf.int = TRUE, conf.level = 0.95, paired = FALSE, exact = FALSE, correct = TRUE)

#Kolmogorov-Smirnov Test:

#Ho: The RI values of Overrepresented anatomical elements in the fossil record and the RI values of Underrepresented anatomical elements in the fossil record are drawn from populations of data with the same cumulative distribution.

#Ha: The RI values of overrepresented anatomical elements in the fossil record and the RI values of underreprensented anatomical elements in the fossil record are drawn from populations of data that have different cumulative distributions. 


ks.test(Overrep, Underrep)

#=======================================================================================

#3. Delta-Statistic (Borges et al., 2019)

#Install packages and set directories
library(ape)
library(phangorn)
library(phytools)
library(devtools)
library(adephylo)
install_github("rgriff23/btw") #BayesTraits Wrapper for R
library("btw")
library(readxl)
library(lattice)

#Running this part of the code requires installing BayesTraits (Meade & Pagel, 2016). For the BayesTraits Wrapper to work in R, the working directory must be set to where you have installed BayesTraits. All your trees and dataframes need to be in this BayesTraits directory as well.
setwd("~/Desktop/Research/Phylogenetic_Signal/BayesTraitsV3.0.2-OSX")

#The first ~200 lines of code are copied directly from Borges et al. (2019) with only one slight modificaiton

#NENTROPY
#returns the node entropies by calculating sum of the state entropies
#prob: matrix of state probabilities

nentropy <- function(prob) {
  
  k              <- ncol(prob)                       #number of states
  prob[prob>1/k] <- prob[prob>1/k]/(1-k) - 1/(1-k)   #state entropies
  tent           <- apply(prob,1,sum)                #node entropy
  
  #correct absolute 0/1
  tent[tent == 0] <- tent[tent == 0] + runif(1,0,1)/10000
  tent[tent == 1] <- tent[tent == 1] - runif(1,0,1)/10000
  
  return(tent)
}

#FUNCTION FOR BAYESIAN INFERENCES
#bayesian inferences on the node entropies 
#l0: rate parameter of the exponential prior distribution
#se: standard deviation of the proposal distribution 
#a:  alpha parameter (beta likelihood)
#b:  beta paramter (beta likelihood)
#x:  node entropies

lpalpha <- function(a,b,x,l0) {          #log posterior alpha
  N  <- length(x)
  lp <- N*(lgamma(a+b)-lgamma(a)) - a*(l0-sum(log(x)))
  return(lp)
}

lpbeta  <- function(a,b,x,l0) {          #log posterior beta
  N  <- length(x)
  lp <- N*(lgamma(a+b)-lgamma(b)) - b*(l0-sum(log(1-x)))
  return(lp)
}

mhalpha <- function(a,b,x,l0,se) {       #metropolis hastings alpha
  a0 <- a
  a1 <- exp(rnorm(1,log(a0),se))
  
  r  <- min(1, exp(lpalpha(a1,b,x,l0) - lpalpha(a0,b,x,l0) ) )
  
  while (is.na(r) == T) {
    a1 <- exp(rnorm(1,log(a0),se))
    r  <- min(1, exp(lpalpha(a1,b,x,l0) - lpalpha(a0,b,x,l0) ) )
  }
  
  if (runif(1) < r) {
    return(a1) 
  } else {
    return(a0)
  }
}

mhbeta  <- function(a,b,x,l0,se) {      #metropolis hastings beta
  b0 <- b
  b1 <- exp(rnorm(1,log(b0),se))
  
  r  <- min(1, exp(lpbeta(a,b1,x,l0) - lpbeta(a,b0,x,l0) ) )
  
  while (is.na(r) == T) {
    b1 <- exp(rnorm(1,log(b0),se))
    r  <- min(1, exp(lpbeta(a,b1,x,l0) - lpbeta(a,b0,x,l0) ) )
  }  
  
  if (runif(1) < r) {
    return(b1)
  } else {
    return(b0)
  }
}

#MCMC
#Markov chain monte carlo scheme using the conditional posteriors of alpha and beta
#alpha: initial value of alpha
#beta: initial values of beta
#x: node entropies
#sim: number of iterations
#thin: controles the number of saved iterations = sim/thin
#burn: number of iterates to burn

emcmc <- function(alpha,beta,x,l0,se,sim,thin,burn) {
  
  usim <- seq(burn,sim,thin)
  gibbs <- matrix(NA,ncol=2,nrow=length(usim))
  p <- 1
  
  for (i in 1:sim) {
    alpha <- mhalpha(alpha,beta,x,l0,se)
    beta  <- mhbeta(alpha,beta,x,l0,se)
    
    if (i == usim[p]) {
      gibbs[p,] <- c(alpha,beta)
      p <- p+1
    }
  }  
  return(gibbs)
}

#RATE MATRIX FOR TRAIT EVOLUTION. K=2 TO 5
ratematrix <- function(pi,rho){
  
  k <- length(pi)
  
  if (k==2){
    r <- c(pi[1]*0     ,pi[2]*rho[1],
           pi[1]*rho[1],pi[2]*0)
  }
  
  if (k==3){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[3],
           pi[1]*rho[2],pi[2]*rho[3],pi[3]*0 )
  }
  
  if (k==4){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],pi[4]*rho[3],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[4],pi[4]*rho[5],
           pi[1]*rho[2],pi[2]*rho[4],pi[3]*0     ,pi[4]*rho[6],
           pi[1]*rho[3],pi[2]*rho[5],pi[3]*rho[6],pi[4]*0 )
  }  
  
  if (k==5){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],pi[4]*rho[3] ,pi[5]*rho[4],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[5],pi[4]*rho[6] ,pi[5]*rho[7],
           pi[1]*rho[2],pi[2]*rho[5],pi[3]*0     ,pi[4]*rho[8] ,pi[5]*rho[9],
           pi[1]*rho[3],pi[2]*rho[6],pi[3]*rho[8],pi[4]*0      ,pi[5]*rho[10],
           pi[1]*rho[4],pi[2]*rho[7],pi[3]*rho[9],pi[4]*rho[10],pi[5]*0)
  }
  
  R <- matrix(r,ncol=k,nrow=k) 
  diag(R) <- -rowSums(R)
  
  return(R)
}

#RTRAIT
#simulates the evolution of a trait in a given tree
# tree: metric-tree
# R: rate matrix
# nstates: number of states

rtrait <- function(tree,R,nstates) {
  
  nspecis <- length(tree$tip.label)
  
  #tree
  edge <- cbind(tree$edge,tree$edge.length)
  
  ancestral <- rep(NA,2*nspecies-1) 
  ancestral[nspecies+1] <- sample(1:nstates,1,prob=pi) 
  
  #rate change
  inode <- nspecies+1
  while (sum(is.na(ancestral)) > 0) {
    
    inode1 <-  edge[which(edge[,1]==inode)[1],2]
    inode2 <-  edge[which(edge[,1]==inode)[2],2]
    bl1 <- edge[which(edge[,1]==inode)[1],3]
    bl2 <- edge[which(edge[,1]==inode)[2],3]
    
    astate <- rep(0,nstates)
    astate[ancestral[inode]] <- 1 
    
    ancestral[inode1] <- sample(1:nstates,1,prob=astate%*%expm(R*bl1))
    ancestral[inode2] <- sample(1:nstates,1,prob=astate%*%expm(R*bl2))
    
    inode <- inode+1
  }
  return(ancestral[1:nspecies])
  
}

delta <- function(tree,lambda0,se,sim,thin,burn) {
  
  #ar below is set equal to the node probabilities inferred from bayestraits. This is the only difference between our analyses and the original script by Borges et al.
  ar <- Node_List
  
  # deletes the complex part whenever it pops up
  if (class(ar[1,1]) == "complex"){
    ar <- Re(ar)
  }
  
  x  <- nentropy(ar)
  mc1    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mc2    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mchain <- rbind(mc1,mc2)
  deltaA <- mean(mchain[,2]/mchain[,1])
  
  return(deltaA)
}




#Import and translate tree
tree<- read.nexus("~/Directory/_Tree.nexus")
write.nexus(tree, file="~/Directory/Tree_Translated.nexus", translate=TRUE)

#tree needs to be re-read into R
tree<- read.nexus("~/Directory/Tree_Translated.nexus")

#drop fossil tips
tree_drop <- drop.tip(tree, tip = c("fossil tips")

#Import matrix. For BayesTraits to read the character matrix, it must be converted to a ".xlsx" file. REMEMBER: any characters that are constant (unchanging character state among all taxa) MUST BE REMOVED. Characters with "?" or "N/A" character states must be changed to "-".
matrix<-read_excel("dataframe_dropped_fossils_and_dropped_characters_with_dash.xlsx", col_names = FALSE)

#Command file for Bayestraits. ##VERY IMPORTANT## The order of these nodes needs to reflect the ascending order of numbered nodes on the tree using nodelabels()
command_vec1 <- c("1", "2", "Iterations 10000000", "Burnin 2000000", "PriorAll uniform 0 1.0", "AddTag Node01_tag Genus_species Genus_species", "AddNode Node01 Node01_tag")

#Assign matrix below. By using these functions, this[1:2] or [c(1,147)], it's possible to run the analysis over different matrices, or poritons of the matrix
Input_Matrix<-matrix

#Assign tree below
Input_Tree<-tree_drop

#Create matrices and lists to hold output 
Deltas<-matrix(nrow=1, ncol=ncol(Input_Matrix)-1)
Nchar<-matrix(nrow=1, ncol=ncol(Input_Matrix)-1)
Node_Values<-list()

#Infer node probabilities and delta statistics for each character
for (l in 1:(ncol(Input_Matrix)-1)){
  Char<-cbind(Input_Matrix[,1],Input_Matrix[,1+l])
  results_1<- bayestraits(Char, Input_Tree, command_vec1)
  results_log<-results_1$Log
  if (ncol(results_log$results)==(Input_Tree$Nnode*2+7)){ # for 2 character states #77
    Values<-results_log$results[,8:ncol(results_log$results)]
    
  } else if (ncol(results_log$results)==(Input_Tree$Nnode*3+12)){ # for 3 character states #117
    Values<-results_log$results[,13:ncol(results_log$results)]
    
  } else if (ncol(results_log$results)==(Input_Tree$Nnode*4+19)){ # for 4 character states #159
    Values<-results_log$results[,20:ncol(results_log$results)]
    
  } else if (ncol(results_log$results)==(Input_Tree$Nnode*5+28)){ # for 5 character states #203
    Values<-results_log$results[,29:ncol(results_log$results)]
  }
  
  
  Mean_Values<-colMeans(Values) #Need to make all numeric for this to work
  if (length(Mean_Values)==(Input_Tree$Nnode*2)) { # for 2 character states #70
    Node_List<-matrix(nrow=length(Mean_Values), ncol=2)
    for (a in seq(from=1, to=length(Mean_Values), by=2)){
      Node_List[a,]<-c(Mean_Values[a], Mean_Values[a+1])
      
    }
  } else if (length(Mean_Values)==(Input_Tree$Nnode*3)){ # for 3 character states #105
    Node_List<-matrix(nrow=length(Mean_Values), ncol=3) 
    for (a in seq(from=1, to=length(Mean_Values), by=3)){
      
      Node_List[a,]<-c(Mean_Values[a], Mean_Values[a+1], Mean_Values[a+2])
      
    }
  }  else if (length(Mean_Values)==(Input_Tree$Nnode*4)){ # for 4 character states #140
    Node_List<-matrix(nrow=length(Mean_Values), ncol=4) 
    for (a in seq(from=1, to=length(Mean_Values), by=4)){
      
      Node_List[a,]<-c(Mean_Values[a], Mean_Values[a+1], Mean_Values[a+2], Mean_Values[a+3])
      
    }
  } else if (length(Mean_Values)==(Input_Tree$Nnode*5)){ # for 5 character states #175
    Node_List<-matrix(nrow=length(Mean_Values), ncol=5) 
    for (a in seq(from=1, to=length(Mean_Values), by=5)){
      
      Node_List[a,]<-c(Mean_Values[a], Mean_Values[a+1], Mean_Values[a+2], Mean_Values[a+3], Mean_Values[a+4])
      
    }
  }
  Node_List<-na.omit(Node_List)
  
  #Below are testing lines to check progress/correct implementation
  print(paste("Character", l)) #Prints Character number
  print(paste(ncol(Node_List),  "character states")) #Prints Number of Character states
  if(sum(Node_List[1,])){ #Should always be TRUE 
    print("TRUE")
  }
  print("x") #Space to make output easier to read
  Nchar[l]<-ncol(Node_List)
  
  #Below are outputs for Node values and Delta statistics for each character
  Node_Values[[l]]<-Node_List
  #row.names(Node_Values[])<-seq(39:73)
  Deltas[l]<-delta(Input_Tree, .1, .5, 10000, 10, 2500)  #100=total number of iterations per character, should be at least 10,000 for final analysis;10= number of samples; 10= burnin (should be closer to 25%)
}


#--------------------------------------------------------------------------------------------------------------------------------------
#GRAPHICAL TOOLS

#Below is optional, but can be used to print nodes on the tree as a pie chart for each character. Only use for a single character.
Plot_Tree<-Input_Tree
plot(Plot_Tree, cex=0.5)
Pie_colors<-c("dodgerblue3", "gold2", "sienna2", "darkorchid3", "aquamarine3")
nodelabels(pie=Node_List, piecol=Pie_colors, cex=0.7, frame="n")
Tip_Labels<-as.data.frame(Plot_Tree$tip.label)
Plot_Tree$tip.label

Char_sorted<-Char[order(match(Char[,1], Tip_Labels[,1])),]

Tip_col<-character(length(Char_sorted))
Tip_col[Char_sorted[,2]=="0"]<-Pie_colors[1]
Tip_col[Char_sorted[,2]=="1"]<-Pie_colors[2]
Tip_col[Char_sorted[,2]=="2"]<-Pie_colors[3]
Tip_col[Char_sorted[,2]=="3"]<-Pie_colors[4]
Tip_col[Char_sorted[,2]=="4"]<-Pie_colors[5]
Tip_col[Char_sorted[,2]=="-"]<-"black"
tiplabels(pch=21, adj=-7, col="black",  bg=Tip_col, cex=1)
}


#Transpose for histogram/other graphical representation of distributions
Delta_values<-t(Deltas)

#Initial boxplot
boxplot(Delta_values, ylim=c(0,6), names="All")

#Example of Anatomical bins for comparison
Tooth_bearing=Delta_values[c(1:6, 10:15, 66:89, 133:178), ]
Cranial=Delta_values[c(7:9, 16:65, 90:132), ]
Axial=Delta_values[179:213, ]
Pectoral=Delta_values[214:230, ]
Pelvic=Delta_values[231:239, ]
Appendicular=Delta_values[240:253, ]
Other_ossifications=Delta_values[254:256, ]
All=Delta_values[1:256, ]

Overrepresented=c(Tooth_bearing, Axial)
Underrepresented=c(Cranial, Pectoral, Pelvic, Appendicular, Other_ossifications)
Total<-c(Sim_Pectoral,Sim_Appendicular)


boxplot(Overrepresented, Underrepresented, ylim=c(0, 6), names=c("Overrepresented","Underrepresented"))

#Plot histograms
histogram(Overrepresented, nint=36, xlim=c(0,0.0), ylim=c(0,35), breaks = do.breaks(endpoints=c(0,100), nint=500000))
histogram(Underrepresented, nint=36, xlim=c(0,0.01), ylim=c(0,35), breaks = do.breaks(endpoints=c(0,100), nint=500000))

#--------------------------------------------------------------------------------------------------------------------------------------

#BOOTSTRAP TESTS FOR DIFFERENCES IN THE MEDIAN VALUES

median(Overrepresented)
median(Underrepresented)

#Calculate difference in medians
Diff_in_real_medians_1<-median(Overrepresented)-median(Underrepresented)

#Bootstrap analysis
Overrepresented_medians<-matrix(c(0),ncol=10000)
Underrepresented_medians<-matrix(c(0),ncol=10000)
Diff_in_medians_1<-matrix(c(0),ncol=10000)
for (m in 1:length(Sim_Appendicular_medians)) {
  Sample_Overrepresented<-sample(Total, size=length(Overrepresented), replace=TRUE)
  Overrepresented_medians[m]<-median(Sample_Overrepresented)
  Sample_Underrepresented<-sample(Total, size=length(Underrepresented), replace=TRUE)
  Underrepresented_medians[m]<-median(Sample_Underrepresented)
  Diff_in_medians_1[m]<-Overrepresented_medians[m]-Underrepresented_medians[m]
}

#Plot up the results
hist(Diff_in_medians_1)
abline(v=Diff_in_real_medians_1)
Values_above<-Diff_in_medians_1[Diff_in_medians_1[1,]>Diff_in_real_medians_1]

#Calculate p-value
p_value<-length(Values_above)/length(Diff_in_medians_1)
p_value


#END
#=======================================================================================




