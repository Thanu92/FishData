#Read the fish dataset in R
fish <- read.table("final_alignment.phylip")
library(tidyverse)
#Rename columns
colnames(fish) <- c("species_name","seq")
#Extract CO1 nucleotides
fishCO1 <- substr(fish$seq,2292,2973)
#Check the class of fishCO1 dataset
class(fishCO1)
#Convert to dataframe
fishCO1 <- data.frame(fishCO1)
#Check the converted dataset
class(fishCO1)
#Bind the column "species name" to the fishCO1 dataframe
fishCO1_with_speciesname <- cbind(fish$species_name,fishCO1)
#Replace "-" with nothing 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)
#Convert to dataframe
df_0 <- as.data.frame(df_0)
#Bind the colum "species name" of fish dataset with new dataframe
df_1 <- cbind(fish$species_name,df_0)
#Remove rows with empty sequences
df_2 <- subset(df_1, df_1$df_0 != "")
#Rename column names
colnames(df_2) <- c("species_name","seq")
#Keep only rows with sequences
fish_multigene <- merge(fish,df_2, by= "species_name")
#New dataframe with first two columns
fish_multigene <- fish_multigene[, 1:2]
#Rename columns
colnames(fish_multigene) <- c("species_name","seq")

#Sample only 100 species for the pilot study
set.seed(1122)
fishsample100 <- fish_multigene[sample(nrow(fish_multigene),100),]
#Check the class
class(fishsample100)
library(ape)
library(phylotools)
#convert datafame to phylip format to use in RAxML
dat2phylip(fishsample100,outfile= "fishsample100.phy")
#convert datafame to fasta format 
library(seqRFLP)
fishsample100_fasta <- dataframe2fas(fishsample100, file = "fishsample100.fasta") 
#dat2fasta(fishsample100, outfile = "fishsample100.fasta")#Not working
#Function to generate different level samples
fishsamples <- function (fishsample_n,y){
  #fishsample_y <- fishsample100[sample(nrow(fishsample100), y),]
  #generate 10 replicatesout without replacement
  fishsample_y <- as.data.frame(replicate(10, fishsample_n[sample(nrow(fishsample_n),y, replace=F),])) 
  return(fishsample_y)
}
set.seed(13)
#Use the function "fishsamples" to generate 80% dataset
fishsample80 <- fishsamples(fishsample100,80)
#Separate dataframes in columns of mass fishsample80 dataframe
replicate_dataframes <- c(fishsample80$V1,fishsample80$V2,fishsample80$V3,fishsample80$V4,fishsample80$V5,fishsample80$V6,fishsample80$V7,fishsample80$V8,fishsample80$V9,fishsample80$V10)
#function to generate dataframes from replicates
df <- function(rep_df){
  t8 <- as.data.frame(rep_df)
  return(t8)
}
# Use the function to generate replicate dataframes 
t9 <- df(replicate_dataframes)
#Separte mass dataframe into a list
lst1 <- lapply(seq(1, ncol(t9), by=2), function(i) 
  t9[i: pmin((i+1), ncol(t9))])
#to get individual dataframes for each replicate
list2env(setNames(lst1, paste0("newdf", seq_along(lst1))),
         envir=.GlobalEnv)
#Function to generate different level samples without replicates and replaces
fishsamples_no_rep <- function (fishsample_m,z){
  fishsample_z <- as.data.frame(fishsample_m[sample(nrow(fishsample_m), z),])
  #fishsample_y <- as.data.frame(replicate(10, fishsample_n[sample(nrow(fishsample_n),y, replace=F),])) 
  return(fishsample_z)
}
#Use the function to genetare phylip files from the dataframe
dff <- function(xx){
  
  return(dat2phylip(xx,outfile = "newphy.phy"))
  
}
#80% 1st dataframe convert to phylip
dff(newdf1)
fishsample60_1 <- fishsamples_no_rep(newdf1,60)
dff(fishsample60_1)
fishsample40_1 <- fishsamples_no_rep(fishsample60_1,40)
dff(fishsample40_1)
fishsample20_1 <- fishsamples_no_rep(fishsample40_1,20)
dff(fishsample20_1)
#80% 2nd dataframe convert to phylip
dff(newdf2)
fishsample60_2 <- fishsamples_no_rep(newdf2,60)
dff(fishsample60_2)
fishsample40_2 <- fishsamples_no_rep(fishsample60_2,40)
dff(fishsample40_2)
fishsample20_2 <- fishsamples_no_rep(fishsample40_2,20)
dff(fishsample20_2)
#80% 3rd dataframe convert to phylip
dff(newdf3)
fishsample60_3 <- fishsamples_no_rep(newdf3,60)
dff(fishsample60_3)
fishsample40_3 <- fishsamples_no_rep(fishsample60_3,40)
dff(fishsample40_3)
fishsample20_3 <- fishsamples_no_rep(fishsample40_3,20)
dff(fishsample20_3)
#80% 4th dataframe convert to phylip
dff(newdf4)
fishsample60_4 <- fishsamples_no_rep(newdf4,60)
dff(fishsample60_4)
fishsample40_4 <- fishsamples_no_rep(fishsample60_4,40)
dff(fishsample40_4)
fishsample20_4 <- fishsamples_no_rep(fishsample40_4,20)
dff(fishsample20_4)
#80% 5th dataframe convert to phylip
dff(newdf5)
fishsample60_5 <- fishsamples_no_rep(newdf5,60)
dff(fishsample60_5)
fishsample40_5 <- fishsamples_no_rep(fishsample60_5,40)
dff(fishsample40_5)
fishsample20_5 <- fishsamples_no_rep(fishsample40_5,20)
dff(fishsample20_5)
#80% 6th dataframe convert to phylip
dff(newdf6)
fishsample60_6 <- fishsamples_no_rep(newdf6,60)
dff(fishsample60_6)
fishsample40_6 <- fishsamples_no_rep(fishsample60_6,40)
dff(fishsample40_6)
fishsample20_6 <- fishsamples_no_rep(fishsample40_6,20)
dff(fishsample20_6)
#80% 7th dataframe convert to phylip
dff(newdf7)
fishsample60_7 <- fishsamples_no_rep(newdf7,60)
dff(fishsample60_7)
fishsample40_7 <- fishsamples_no_rep(fishsample60_7,40)
dff(fishsample40_7)
fishsample20_7 <- fishsamples_no_rep(fishsample40_7,20)
dff(fishsample20_7)
#80% 8th dataframe convert to phylip
dff(newdf8)
fishsample60_8 <- fishsamples_no_rep(newdf8,60)
dff(fishsample60_8)
fishsample40_8 <- fishsamples_no_rep(fishsample60_8,40)
dff(fishsample40_8)
fishsample20_8 <- fishsamples_no_rep(fishsample40_8,20)
dff(fishsample20_8)
#80% 9th dataframe convert to phylip
dff(newdf9)
fishsample60_9 <- fishsamples_no_rep(newdf9,60)
dff(fishsample60_9)
fishsample40_9 <- fishsamples_no_rep(fishsample60_9,40)
dff(fishsample40_9)
fishsample20_9 <- fishsamples_no_rep(fishsample40_9,20)
dff(fishsample20_9)
#80% 10th dataframe convert to phylip
dff(newdf10)
fishsample60_10 <- fishsamples_no_rep(newdf10,60)
dff(fishsample60_10)
fishsample40_10 <- fishsamples_no_rep(fishsample60_10,40)
dff(fishsample40_10)
fishsample20_10 <- fishsamples_no_rep(fishsample40_10,20)
dff(fishsample20_10)
#function of finding missing data
missing_data <- function(level_df){
  data_missed <- data.frame(dplyr::setdiff(fishsample100,level_df))
  return(data_missed)
}
#Assign column names
colnames(newdf2) <- c("species_name","seq")
colnames(newdf3) <- c("species_name","seq")
colnames(newdf4) <- c("species_name","seq")
colnames(newdf5) <- c("species_name","seq")
colnames(newdf6) <- c("species_name","seq")
colnames(newdf7) <- c("species_name","seq")
colnames(newdf8) <- c("species_name","seq")
colnames(newdf9) <- c("species_name","seq")
colnames(newdf10) <- c("species_name","seq")
colnames(fishsample60_2) <- c("species_name","seq")
colnames(fishsample60_3) <- c("species_name","seq")
colnames(fishsample60_4) <- c("species_name","seq")
colnames(fishsample60_5) <- c("species_name","seq")
colnames(fishsample60_6) <- c("species_name","seq")
colnames(fishsample60_7) <- c("species_name","seq")
colnames(fishsample60_8) <- c("species_name","seq")
colnames(fishsample60_9) <- c("species_name","seq")
colnames(fishsample60_10) <- c("species_name","seq")
colnames(fishsample40_2) <- c("species_name","seq")
colnames(fishsample40_3) <- c("species_name","seq")
colnames(fishsample40_4) <- c("species_name","seq")
colnames(fishsample40_5) <- c("species_name","seq")
colnames(fishsample40_6) <- c("species_name","seq")
colnames(fishsample40_7) <- c("species_name","seq")
colnames(fishsample40_8) <- c("species_name","seq")
colnames(fishsample40_9) <- c("species_name","seq")
colnames(fishsample40_10) <- c("species_name","seq")
colnames(fishsample20_2) <- c("species_name","seq")
colnames(fishsample20_3) <- c("species_name","seq")
colnames(fishsample20_4) <- c("species_name","seq")
colnames(fishsample20_5) <- c("species_name","seq")
colnames(fishsample20_6) <- c("species_name","seq")
colnames(fishsample20_7) <- c("species_name","seq")
colnames(fishsample20_8) <- c("species_name","seq")
colnames(fishsample20_9) <- c("species_name","seq")
colnames(fishsample20_10) <- c("species_name","seq")

#Use the function of finding missing data
miss_80_1 <- missing_data(newdf1)
miss_80_2 <- missing_data(newdf2)
miss_80_3 <- missing_data(newdf3)
miss_80_4 <- missing_data(newdf4)
miss_80_5 <- missing_data(newdf5)
miss_80_6 <- missing_data(newdf6)
miss_80_7 <- missing_data(newdf7)
miss_80_8 <- missing_data(newdf8)
miss_80_9 <- missing_data(newdf9)
miss_80_10 <- missing_data(newdf10)
miss_60_1 <- missing_data(fishsample60_1)
miss_60_2 <- missing_data(fishsample60_2)
miss_60_3 <- missing_data(fishsample60_3)
miss_60_4 <- missing_data(fishsample60_4)
miss_60_5 <- missing_data(fishsample60_5)
miss_60_6 <- missing_data(fishsample60_6)
miss_60_7 <- missing_data(fishsample60_7)
miss_60_8 <- missing_data(fishsample60_8)
miss_60_9 <- missing_data(fishsample60_9)
miss_60_10 <- missing_data(fishsample60_10)
miss_40_1 <- missing_data(fishsample40_1)
miss_40_2 <- missing_data(fishsample40_2)
miss_40_3 <- missing_data(fishsample40_3)
miss_40_4 <- missing_data(fishsample40_4)
miss_40_5 <- missing_data(fishsample40_5)
miss_40_6 <- missing_data(fishsample40_6)
miss_40_7 <- missing_data(fishsample40_7)
miss_40_8 <- missing_data(fishsample40_8)
miss_40_9 <- missing_data(fishsample40_9)
miss_40_10 <- missing_data(fishsample40_10)
miss_20_1 <- missing_data(fishsample20_1)
miss_20_2 <- missing_data(fishsample20_2)
miss_20_3 <- missing_data(fishsample20_3)
miss_20_4 <- missing_data(fishsample20_4)
miss_20_5 <- missing_data(fishsample20_5)
miss_20_6 <- missing_data(fishsample20_6)
miss_20_7 <- missing_data(fishsample20_7)
miss_20_8 <- missing_data(fishsample20_8)
miss_20_9 <- missing_data(fishsample20_9)
miss_20_10 <- missing_data(fishsample20_10)

#Extract the range of COI gene
miss_80_1_co1 <- substr(miss_80_1$seq,2292,2973)
#Check the class
class(miss_80_1_co1)
#transfered to a dataframe
miss_80_1_co1 <- data.frame(miss_80_1_co1)
#Check the class
class(miss_80_1_co1)
#Bind the column of species name
miss_80_1_co1_with_speciesname <- cbind(miss_80_1$species_name,miss_80_1_co1)
#Assign column names
colnames(miss_80_1_co1_with_speciesname) <- c("species_name","seq")
#Transfered to the phylip format
dat2phylip(miss_80_1_co1_with_speciesname,outfile= "miss_80_1_co1.phy")
#Combins rows
miss_rbind_80_1 <- rbind(miss_80_1_co1_with_speciesname,newdf1)
#class(miss_rbind_80_1)
dat2phylip(miss_rbind_80_1,outfile= "rbind_80miss20bb_1.phy")
fishsamplerbind_80miss20bb_1_fasta <- dataframe2fas(miss_rbind_80_1,file= "rbind_80miss20bb_1.fasta")

#sum(is.na(fishCO1_miss_80_1_speciesname$fishCO1_miss_80_1))
#miss_rbind_80_1$seq <- DNAStringSet((miss_rbind_80_1$seq)) #covnert the class of the sequences into Biostrings
#class(miss_rbind_80_1$seq) #check the class

#ab <- DNAStringSet(muscle::muscle(miss_rbind_80_1$seq, log = "log.tx", verbose = T))

#Extract the range of COI gene
fishsample100_co1 <- substr(fishsample100$seq,2292,2973)
#Check the class
class(fishsample100_co1)
#Convert to dataframe
fishsample100_co1 <- data.frame(fishsample100_co1)
#Check the class
class(fishsample100_co1)
#Combine columns of COI dataframe and fishsample100 dataframe to get the sepecies name
fishsample100_co1_with_speciesname <- cbind(fishsample100$species_name,fishsample100_co1)
#Assign column names
colnames(miss_80_1_co1_with_speciesname) <- c("species_name","seq")
#Convert to phylip format
dat2phylip(fishsample100_co1_with_speciesname,outfile= "fishsample100_co1_with_speciesname.phy")
fishsample100_co1_fasta <- dataframe2fas(fishsample100_co1_with_speciesname, file = "fishsample100_co1.fasta") 
#text_file <- read.csv("twoTrees")

#Generate Complete tree and different level backbone trees using RAxML
#Visualize generated backbone trees using FigTree
#Phylogenetic placement using RAxMLEPA
#Open jplace files generated in RAxMLEPA using iTol online tool
#Change the format of jlace into newick format usig iToL online tool
#Export the newick format


library(phangorn)
#Read trees into R generated from RAxML
a <- read.tree("RAxML_parsimonyTree.complete100")
b <- read.tree("20_1_new")   
c <- read.tree("20_2_new") 
d <- read.tree("20_3_new") 
e <- read.tree("20_4_new") 
f <- read.tree("20_5_new") 
g <- read.tree("20_6_new")   
h <- read.tree("20_7_new")     
i <- read.tree("20_8_new") 
j <- read.tree("20_9_new") 
k <- read.tree("20_10_new") 
l <- read.tree("40_1_new") 
m <- read.tree("40_2_new") 
n <- read.tree("40_3_new") 
o <- read.tree("40_4_new") 
p <- read.tree("40_5_new") 
q <- read.tree("40_6_new") 
r <- read.tree("40_7_new") 
s <- read.tree("40_8_new") 
t <- read.tree("40_9_new") 
u <- read.tree("40_10_new") 
v <- read.tree("60_1_new") 
w <- read.tree("60_2_new") 
x <- read.tree("60_3_new") 
y <- read.tree("60_4_new") 
z <- read.tree("60_5_new") 
aa <- read.tree("60_6_new") 
bb <- read.tree("60_7_new") 
cc <- read.tree("60_8_new")
dd <- read.tree("60_9_new")
ee <- read.tree("60_10_new")
ff <- read.tree("80_1_new")
gg <- read.tree("80_2_new")
hh <- read.tree("80_3_new")
ii <- read.tree("80_4_new")
jj <- read.tree("80_5_new")
kk <- read.tree("80_6_new")
ll <- read.tree("80_7_new")
mm <- read.tree("80_8_new")
nn <- read.tree("80_9_new")
oo <- read.tree("80_10_new")
#Robinson Fould distance (Depends on topology of trees)
RF_20_1 <- RF.dist(a,b,normalize = TRUE)
#path distance (Path difference)
P_20_1  <- path.dist(a,b)
RF_20_2 <- RF.dist(a,c,normalize = TRUE)
P_20_2 <- path.dist(a,c)
RF_20_3 <- RF.dist(a,d,normalize = TRUE)
P_20_3 <- path.dist(a,d)
RF_20_4 <- RF.dist(a,e,normalize = TRUE)
P_20_4 <- path.dist(a,e)
RF_20_5 <- RF.dist(a,f,normalize = TRUE)
P_20_5 <- path.dist(a,f)
RF_20_6 <- RF.dist(a,g,normalize = TRUE)
P_20_6 <- path.dist(a,g)
RF_20_7 <- RF.dist(a,h,normalize = TRUE)
P_20_7 <- path.dist(a,h)
RF_20_8 <- RF.dist(a,i,normalize = TRUE)
P_20_8 <- path.dist(a,i)
RF_20_9 <- RF.dist(a,j,normalize = TRUE)
P_20_9 <- path.dist(a,j)
RF_20_10 <- RF.dist(a,k,normalize = TRUE)
P_20_10 <- path.dist(a,k)
RF_40_1 <- RF.dist(a,l,normalize = TRUE)
P_40_1 <- path.dist(a,l)
RF_40_2 <- RF.dist(a,m,normalize = TRUE)
P_40_2 <- path.dist(a,m)
RF_40_3 <- RF.dist(a,n,normalize = TRUE)
P_40_3 <- path.dist(a,n)
RF_40_4 <- RF.dist(a,o,normalize = TRUE)
P_40_4 <- path.dist(a,o)
RF_40_5 <- RF.dist(a,p,normalize = TRUE)
P_40_5 <- path.dist(a,p)
RF_40_6 <- RF.dist(a,q,normalize = TRUE)
P_40_6 <- path.dist(a,q)
RF_40_7 <- RF.dist(a,r,normalize = TRUE)
P_40_7 <- path.dist(a,r)
RF_40_8 <- RF.dist(a,s,normalize = TRUE)
P_40_8 <- path.dist(a,s)
RF_40_9 <- RF.dist(a,t,normalize = TRUE)
P_40_9 <- path.dist(a,t)
RF_40_10 <- RF.dist(a,u,normalize = TRUE)
P_40_10 <- path.dist(a,u)
RF_60_1 <- RF.dist(a,v,normalize = TRUE)
P_60_1 <- path.dist(a,v)
RF_60_2 <- RF.dist(a,w,normalize = TRUE)
P_60_2 <- path.dist(a,w)
RF_60_3 <- RF.dist(a,x,normalize = TRUE)
P_60_3 <- path.dist(a,x)
RF_60_4 <- RF.dist(a,y,normalize = TRUE)
P_60_4 <- path.dist(a,y)
RF_60_5 <- RF.dist(a,z,normalize = TRUE)
P_60_5 <- path.dist(a,z)
RF_60_6 <- RF.dist(a,aa,normalize = TRUE)
P_60_6 <- path.dist(a,aa)
RF_60_7 <- RF.dist(a,bb,normalize = TRUE)
P_60_7 <- path.dist(a,bb)
RF_60_8 <- RF.dist(a,cc,normalize = TRUE)
P_60_8 <- path.dist(a,cc)
RF_60_9 <- RF.dist(a,dd,normalize = TRUE)
P_60_9 <- path.dist(a,dd)
RF_60_10 <- RF.dist(a,ee,normalize = TRUE)
P_60_10 <- path.dist(a,ee)
RF_80_1 <- RF.dist(a,ff,normalize = TRUE)
P_80_1 <- path.dist(a,ff)
RF_80_2 <- RF.dist(a,gg,normalize = TRUE)
P_80_2 <- path.dist(a,gg)
RF_80_3 <- RF.dist(a,hh,normalize = TRUE)
P_80_3 <- path.dist(a,hh)
RF_80_4 <- RF.dist(a,ii,normalize = TRUE)
P_80_4 <- path.dist(a,ii)
RF_80_5 <- RF.dist(a,jj,normalize = TRUE)
P_80_5 <- path.dist(a,jj)
RF_80_6 <- RF.dist(a,kk,normalize = TRUE)
P_80_6 <- path.dist(a,kk)
RF_80_7 <- RF.dist(a,ll,normalize = TRUE)
P_80_7 <- path.dist(a,ll)
RF_80_8 <- RF.dist(a,mm,normalize = TRUE)
P_80_8 <- path.dist(a,mm)
RF_80_9 <- RF.dist(a,nn,normalize = TRUE)
P_80_9 <- path.dist(a,nn)
RF_80_10 <- RF.dist(a,oo,normalize = TRUE)
P_80_10 <- path.dist(a,oo)
#plot Robinson Fould distance
GRF <- c(20, 40, 60, 80)
GRF_1 <- c(RF_20_1, RF_40_1, RF_60_1, RF_80_1)
GRF_2 <- c(RF_20_2, RF_40_2, RF_60_2, RF_80_2)
GRF_3 <- c(RF_20_3, RF_40_3, RF_60_3, RF_80_3)
GRF_4 <- c(RF_20_4, RF_40_4, RF_60_4, RF_80_4)
GRF_5 <-  c(RF_20_5, RF_40_5, RF_60_5, RF_80_5)
GRF_6 <- c(RF_20_6, RF_40_6, RF_60_6, RF_80_6)
GRF_7 <- c(RF_20_7, RF_40_7, RF_60_7, RF_80_7)
GRF_8 <- c(RF_20_8, RF_40_8, RF_60_8, RF_80_8)
GRF_9 <- c(RF_20_9, RF_40_9, RF_60_9, RF_80_9)
GRF_10 <- c(RF_20_10, RF_40_10, RF_60_10, RF_80_10)

GP_1 <- c(P_20_1, P_40_1, P_60_1, P_80_1)
GP_2 <- c(P_20_2, P_40_2, P_60_2, P_80_2)
GP_3 <- c(P_20_3, P_40_3, P_60_3, P_80_3)
GP_4 <- c(P_20_4, P_40_4, P_60_4, P_80_4)
GP_5 <-  c(P_20_5, P_40_5, P_60_5, P_80_5)
GP_6 <- c(P_20_6, P_40_6, P_60_6, P_80_6)
GP_7 <- c(P_20_7, P_40_7, P_60_7, P_80_7)
GP_8 <- c(P_20_8, P_40_8, P_60_8, P_80_8)
GP_9 <- c(P_20_9, P_40_9, P_60_9, P_80_9)
GP_10 <- c(P_20_10, P_40_10, P_60_10, P_80_10)


dffff <- data.frame(GRF,GRF_1,GRF_2,GRF_3,GRF_4,GRF_5,GRF_6,GRF_7,GRF_8,GRF_9,GRF_10)
dffffF <- data.frame(GRF,GP_1,GP_2,GP_3,GP_4,GP_5,GP_6,GP_7,GP_8,GP_9,GP_10)
ggplot(dffff  , aes(GRF)) +
  scale_y_reverse() + xlab("The level of backbone(%)")+ylab("RF Distance")+labs(title="The plot of RF distance vs. the level of backbone")+
  geom_line(aes(y=GRF_1),
            colour="red") +
  geom_line(aes(y=GRF_2),
            colour="green") +
  geom_line(aes(y=GRF_3),
            colour="blue") +
  geom_line(aes(y=GRF_4),
            colour="yellow") +
  geom_line(aes(y=GRF_5),
            colour="orange") +
  geom_line(aes(y=GRF_6),
            colour="black")+
  geom_line(aes(y=GRF_7),
            colour="purple") +
  geom_line(aes(y=GRF_8),
            colour="pink") +
  geom_line(aes(y=GRF_9),
            colour="brown") +
  geom_line(aes(y=GRF_10),
            colour="violet")
#plot path distances
ggplot(dffffF  , aes(GRF)) +
  scale_y_reverse() + xlab("The level of backbone(%)")+ylab("Path Distance")+labs(title="The plot of Path distance vs. the level of backbone")+
  geom_line(aes(y=GP_1),
            colour="red") +
  geom_line(aes(y=GP_2),
            colour="green") +
  geom_line(aes(y=GP_3),
            colour="blue") +
  geom_line(aes(y=GP_4),
            colour="yellow") +
  geom_line(aes(y=GP_5),
            colour="orange") +
  geom_line(aes(y=GP_6),
            colour="black")+
  geom_line(aes(y=GP_7),
            colour="purple") +
  geom_line(aes(y=GP_8),
            colour="pink") +
  geom_line(aes(y=GP_9),
            colour="brown") +
  geom_line(aes(y=GP_10),
            colour="violet")

library(ggplot2)
library(dplyr)

library(ape)
#Newick to distancematrix
DM_100 <- cophenetic(a)/ max(cophenetic(a))
DM_20_1 <- cophenetic(b)/ max(cophenetic(b))
DM_20_2 <- cophenetic(c)/ max(cophenetic(c))
DM_20_3 <- cophenetic(d)/ max(cophenetic(d))
DM_20_4 <- cophenetic(e)/ max(cophenetic(e))
DM_20_5 <- cophenetic(f)/ max(cophenetic(f))
DM_20_6 <- cophenetic(g)/ max(cophenetic(g))
DM_20_7 <- cophenetic(h)/ max(cophenetic(h))
DM_20_8 <- cophenetic(i)/ max(cophenetic(i))
DM_20_9 <- cophenetic(j)/ max(cophenetic(j))
DM_20_10 <- cophenetic(k)/ max(cophenetic(k))
DM_40_1 <- cophenetic(l)/ max(cophenetic(l))
DM_40_2 <- cophenetic(m)/ max(cophenetic(m))
DM_40_3 <- cophenetic(n)/ max(cophenetic(n))
DM_40_4 <- cophenetic(o)/ max(cophenetic(o))
DM_40_5 <- cophenetic(p)/ max(cophenetic(p))
DM_40_6 <- cophenetic(q)/ max(cophenetic(q))
DM_40_7 <- cophenetic(r)/ max(cophenetic(r))
DM_40_8 <- cophenetic(s)/ max(cophenetic(s))
DM_40_9 <- cophenetic(t)/ max(cophenetic(t))
DM_40_10 <- cophenetic(u)/ max(cophenetic(u))
DM_60_1 <- cophenetic(v)/ max(cophenetic(v))
DM_60_2 <- cophenetic(w)/ max(cophenetic(w))
DM_60_3 <- cophenetic(x)/ max(cophenetic(x))
DM_60_4 <- cophenetic(y)/ max(cophenetic(y))
DM_60_5 <- cophenetic(z)/ max(cophenetic(z))
DM_60_6 <- cophenetic(aa)/ max(cophenetic(aa))
DM_60_7 <- cophenetic(bb)/ max(cophenetic(bb))
DM_60_8 <- cophenetic(cc)/ max(cophenetic(cc))
DM_60_9 <- cophenetic(dd)/ max(cophenetic(dd))
DM_60_10 <- cophenetic(ee)/ max(cophenetic(ee))
DM_80_1 <- cophenetic(ff)/ max(cophenetic(ff))
DM_80_2 <- cophenetic(gg)/ max(cophenetic(gg))
DM_80_3 <- cophenetic(hh)/ max(cophenetic(hh))
DM_80_4 <- cophenetic(ii)/ max(cophenetic(ii))
DM_80_5 <- cophenetic(jj)/ max(cophenetic(jj))
DM_80_6 <- cophenetic(kk)/ max(cophenetic(kk))
DM_80_7 <- cophenetic(ll)/ max(cophenetic(ll))
DM_80_8 <- cophenetic(mm)/ max(cophenetic(mm))
DM_80_9 <- cophenetic(nn)/ max(cophenetic(nn))
DM_80_10 <- cophenetic(oo)/ max(cophenetic(oo))
#function to generate CADM results
CADM <- function(DM){
  #add two matrics together
  AB <- rbind(DM_100, DM)
  cadm_global <- CADM.global(AB, 2, 100)$congruence_analysis[1,]
  return(cadm_global)
}
#cadm_AB <- CADM(DM_B)
#cadm_AB$congruence_analysis
#Use CADM function to compare backbone trees with the complete tree
CADM(DM_20_1)
CADM(DM_20_2)
CADM(DM_20_3)
CADM(DM_20_4)
CADM(DM_20_5)
CADM(DM_20_6)
CADM(DM_20_7)
CADM(DM_20_8)
CADM(DM_20_9)
CADM(DM_20_10)
CADM(DM_40_1)
CADM(DM_40_2)
CADM(DM_40_3)
CADM(DM_40_4)
CADM(DM_40_5)
CADM(DM_40_6)
CADM(DM_40_7)
CADM(DM_40_8)
CADM(DM_40_9)
CADM(DM_40_10)
CADM(DM_60_1)
CADM(DM_60_2)
CADM(DM_60_3)
CADM(DM_60_4)
CADM(DM_60_5)
CADM(DM_60_6)
CADM(DM_60_7)
CADM(DM_60_8)
CADM(DM_60_9)
CADM(DM_60_10)
CADM(DM_80_1)
CADM(DM_80_2)
CADM(DM_80_3)
CADM(DM_80_4)
CADM(DM_80_5)
CADM(DM_80_6)
CADM(DM_80_7)
CADM(DM_80_8)
CADM(DM_80_9)
CADM(DM_80_10)
#plot congruence among distance matric
CADM_1 <- c(CADM(DM_20_1), CADM(DM_40_1), CADM(DM_60_1), CADM(DM_80_1))
CADM_2 <- c(CADM(DM_20_2), CADM(DM_40_2), CADM(DM_60_2), CADM(DM_80_2))
CADM_3 <- c(CADM(DM_20_3), CADM(DM_40_3), CADM(DM_60_3), CADM(DM_80_3))
CADM_4 <- c(CADM(DM_20_4), CADM(DM_40_4), CADM(DM_60_4), CADM(DM_80_4))
CADM_5 <- c(CADM(DM_20_5), CADM(DM_40_5), CADM(DM_60_5), CADM(DM_80_5))
CADM_6 <- c(CADM(DM_20_6), CADM(DM_40_6), CADM(DM_60_6), CADM(DM_80_6))
CADM_7 <- c(CADM(DM_20_7), CADM(DM_40_7), CADM(DM_60_7), CADM(DM_80_7))
CADM_8 <- c(CADM(DM_20_8), CADM(DM_40_8), CADM(DM_60_8), CADM(DM_80_8))
CADM_9 <- c(CADM(DM_20_9), CADM(DM_40_9), CADM(DM_60_9), CADM(DM_80_9))
CADM_10 <- c(CADM(DM_20_10), CADM(DM_40_10), CADM(DM_60_10), CADM(DM_80_10))
df_CADM <- data.frame(GRF,CADM_1,CADM_2,CADM_3,CADM_4,CADM_5,CADM_6,CADM_7,CADM_8,CADM_9,CADM_10)
ggplot(df_CADM  , aes(GRF)) + xlim(20,80) + ylim(0,1) +
  xlab("The level of backbone(%)")+ylab("Congruence Analysis")+labs(title="The plot of Congruence Analysis vs. the level of backbone")+
  geom_line(aes(y=CADM_1),
            colour="red") +
  geom_line(aes(y=CADM_2),
            colour="green") +
  geom_line(aes(y=CADM_3),
            colour="blue") +
  geom_line(aes(y=CADM_4),
            colour="yellow") +
  geom_line(aes(y=CADM_5),
            colour="orange") +
  geom_line(aes(y=CADM_6),
            colour="black")+
  geom_line(aes(y=CADM_7),
            colour="purple") +
  geom_line(aes(y=CADM_8),
            colour="pink") +
  geom_line(aes(y=CADM_9),
            colour="brown") +
  geom_line(aes(y=CADM_10),
            colour="violet")

#install.packages('TreeDist')
library(TreeDist)
D_20_1 <- TreeDistance(a,b)
D_20_2 <- TreeDistance(a,c)
D_20_3 <- TreeDistance(a,d)
D_20_4 <- TreeDistance(a,e)
D_20_5 <- TreeDistance(a,f)
D_20_6 <- TreeDistance(a,g)
D_20_7 <- TreeDistance(a,h)
D_20_8 <- TreeDistance(a,i)
D_20_9 <- TreeDistance(a,j)
D_20_10 <- TreeDistance(a,k)
D_40_1 <- TreeDistance(a,l)
D_40_2 <- TreeDistance(a,m)
D_40_3 <- TreeDistance(a,n)
D_40_4 <- TreeDistance(a,o)
D_40_5 <- TreeDistance(a,p)
D_40_6 <- TreeDistance(a,q)
D_40_7 <- TreeDistance(a,r)
D_40_8 <- TreeDistance(a,s)
D_40_9 <- TreeDistance(a,t)
D_40_10 <- TreeDistance(a,u)
D_60_1 <- TreeDistance(a,v)
D_60_2 <- TreeDistance(a,w)
D_60_3 <- TreeDistance(a,x)
D_60_4 <- TreeDistance(a,y)
D_60_5 <- TreeDistance(a,z)
D_60_6 <- TreeDistance(a,aa)
D_60_7 <- TreeDistance(a,bb)
D_60_8 <- TreeDistance(a,cc)
D_60_9 <- TreeDistance(a,dd)
D_60_10 <- TreeDistance(a,ee)
D_80_1 <- TreeDistance(a,ff)
D_80_2 <- TreeDistance(a,gg)
D_80_3 <- TreeDistance(a,hh)
D_80_4 <- TreeDistance(a,ii)
D_80_5 <- TreeDistance(a,jj)
D_80_6 <- TreeDistance(a,kk)
D_80_7 <- TreeDistance(a,ll)
D_80_8 <- TreeDistance(a,mm)
D_80_9 <- TreeDistance(a,nn)
D_80_10 <- TreeDistance(a,oo)
#Plot treedistances
GD_1 <- c(D_20_1, D_40_1, D_60_1, D_80_1)
GD_2 <- c(D_20_2, D_40_2, D_60_2, D_80_2)
GD_3 <- c(D_20_3, D_40_3, D_60_3, D_80_3)
GD_4 <- c(D_20_4, D_40_4, D_60_4, D_80_4)
GD_5 <-  c(D_20_5, D_40_5, D_60_5, D_80_5)
GD_6 <- c(D_20_6, D_40_6, D_60_6, D_80_6)
GD_7 <- c(D_20_7, D_40_7, D_60_7, D_80_7)
GD_8 <- c(D_20_8, D_40_8, D_60_8, D_80_8)
GD_9 <- c(D_20_9, D_40_9, D_60_9, D_80_9)
GD_10 <- c(D_20_10, D_40_10, D_60_10, D_80_10)

df_GD <- data.frame(GRF,GD_1,GD_2,GD_3,GD_4,GD_5,GD_6,GD_7,GD_8,GD_9,GD_10)
ggplot(df_GD  , aes(GRF)) +
  scale_y_reverse() + xlab("The level of backbone(%)")+ylab("Tree Distance")+labs(title="The plot of Tree distance vs. the level of backbone")+
  geom_line(aes(y=GD_1),
            colour="red") +
  geom_line(aes(y=GD_2),
            colour="green") +
  geom_line(aes(y=GD_3),
            colour="blue") +
  geom_line(aes(y=GD_4),
            colour="yellow") +
  geom_line(aes(y=GD_5),
            colour="orange") +
  geom_line(aes(y=GD_6),
            colour="black")+
  geom_line(aes(y=GD_7),
            colour="purple") +
  geom_line(aes(y=GD_8),
            colour="pink") +
  geom_line(aes(y=GD_9),
            colour="brown") +
  geom_line(aes(y=GD_10),
            colour="violet")
#Calculate SharedPhylogeneticInfo
SPI_20_1 <- SharedPhylogeneticInfo(a, b,normalize = TRUE)
SPI_20_2 <- SharedPhylogeneticInfo(a, c,normalize = TRUE)
SPI_20_3 <- SharedPhylogeneticInfo(a, d,normalize = TRUE)
SPI_20_4 <- SharedPhylogeneticInfo(a, e,normalize = TRUE)
SPI_20_5 <- SharedPhylogeneticInfo(a, f,normalize = TRUE)
SPI_20_6 <- SharedPhylogeneticInfo(a, g,normalize = TRUE)
SPI_20_7 <- SharedPhylogeneticInfo(a, h,normalize = TRUE)
SPI_20_8 <- SharedPhylogeneticInfo(a, i,normalize = TRUE)
SPI_20_9 <- SharedPhylogeneticInfo(a, j,normalize = TRUE)
SPI_20_10 <- SharedPhylogeneticInfo(a, k,normalize = TRUE)
SPI_40_1 <- SharedPhylogeneticInfo(a, l,normalize = TRUE)
SPI_40_2 <- SharedPhylogeneticInfo(a, m,normalize = TRUE)
SPI_40_3 <- SharedPhylogeneticInfo(a, n,normalize = TRUE)
SPI_40_4 <- SharedPhylogeneticInfo(a, o,normalize = TRUE)
SPI_40_5 <- SharedPhylogeneticInfo(a, p,normalize = TRUE)
SPI_40_6 <- SharedPhylogeneticInfo(a, q,normalize = TRUE)
SPI_40_7 <- SharedPhylogeneticInfo(a, r,normalize = TRUE)
SPI_40_8 <- SharedPhylogeneticInfo(a, s,normalize = TRUE)
SPI_40_9 <- SharedPhylogeneticInfo(a, t,normalize = TRUE)
SPI_40_10 <- SharedPhylogeneticInfo(a, u,normalize = TRUE)
SPI_60_1 <- SharedPhylogeneticInfo(a, v,normalize = TRUE)
SPI_60_2 <- SharedPhylogeneticInfo(a, w,normalize = TRUE)
SPI_60_3 <- SharedPhylogeneticInfo(a, x,normalize = TRUE)
SPI_60_4 <- SharedPhylogeneticInfo(a, y,normalize = TRUE)
SPI_60_5 <- SharedPhylogeneticInfo(a, z,normalize = TRUE)
SPI_60_6 <- SharedPhylogeneticInfo(a, aa,normalize = TRUE)
SPI_60_7 <- SharedPhylogeneticInfo(a, bb,normalize = TRUE)
SPI_60_8 <- SharedPhylogeneticInfo(a, cc,normalize = TRUE)
SPI_60_9 <- SharedPhylogeneticInfo(a, dd,normalize = TRUE)
SPI_60_10 <- SharedPhylogeneticInfo(a, ee,normalize = TRUE)
SPI_80_1 <- SharedPhylogeneticInfo(a, ff,normalize = TRUE)
SPI_80_2 <- SharedPhylogeneticInfo(a, gg,normalize = TRUE)
SPI_80_3 <- SharedPhylogeneticInfo(a, hh,normalize = TRUE)
SPI_80_4 <- SharedPhylogeneticInfo(a, ii,normalize = TRUE)
SPI_80_5 <- SharedPhylogeneticInfo(a, jj,normalize = TRUE)
SPI_80_6 <- SharedPhylogeneticInfo(a, kk,normalize = TRUE)
SPI_80_7 <- SharedPhylogeneticInfo(a, ll,normalize = TRUE)
SPI_80_8 <- SharedPhylogeneticInfo(a, mm,normalize = TRUE)
SPI_80_9 <- SharedPhylogeneticInfo(a, nn,normalize = TRUE)
SPI_80_10 <- SharedPhylogeneticInfo(a, oo,normalize = TRUE)
#Plot SharedPhylogeneticInfo
GSPI_1 <- c(SPI_20_1, SPI_40_1, SPI_60_1, SPI_80_1)
GSPI_2 <- c(SPI_20_2, SPI_40_2, SPI_60_2, SPI_80_2)
GSPI_3 <- c(SPI_20_3, SPI_40_3, SPI_60_3, SPI_80_3)
GSPI_4 <- c(SPI_20_4, SPI_40_4, SPI_60_4, SPI_80_4)
GSPI_5 <-  c(SPI_20_5, SPI_40_5, SPI_60_5, SPI_80_5)
GSPI_6 <- c(SPI_20_6, SPI_40_6, SPI_60_6, SPI_80_6)
GSPI_7 <- c(SPI_20_7, SPI_40_7, SPI_60_7, SPI_80_7)
GSPI_8 <- c(SPI_20_8, SPI_40_8, SPI_60_8, SPI_80_8)
GSPI_9 <- c(SPI_20_9, SPI_40_9, SPI_60_9, SPI_80_9)
GSPI_10 <- c(SPI_20_10, SPI_40_10, SPI_60_10, SPI_80_10)

df_GSPI <- data.frame(GRF,GSPI_1,GSPI_2,GSPI_3,GSPI_4,GSPI_5,GSPI_6,GSPI_7,GSPI_8,GSPI_9,GSPI_10)
ggplot(df_GSPI  , aes(GRF)) +
   xlab("The level of backbone(%)")+ylab("Shared Phylogenetic Infomation")+labs(title="The plot of Shared Phylogenetic Infomation vs. the level of backbone")+
  geom_line(aes(y=GSPI_1),
            colour="red") +
  geom_line(aes(y=GSPI_2),
            colour="green") +
  geom_line(aes(y=GSPI_3),
            colour="blue") +
  geom_line(aes(y=GSPI_4),
            colour="yellow") +
  geom_line(aes(y=GSPI_5),
            colour="orange") +
  geom_line(aes(y=GSPI_6),
            colour="black")+
  geom_line(aes(y=GSPI_7),
            colour="purple") +
  geom_line(aes(y=GSPI_8),
            colour="pink") +
  geom_line(aes(y=GSPI_9),
            colour="brown") +
  geom_line(aes(y=GSPI_10),
            colour="violet")
#Calculate MutualClusteringInfo
MCI_20_1 <- MutualClusteringInfo(a, b,normalize = TRUE)
MCI_20_2 <- MutualClusteringInfo(a, c,normalize = TRUE)
MCI_20_3 <- MutualClusteringInfo(a, d,normalize = TRUE)
MCI_20_4 <- MutualClusteringInfo(a, e,normalize = TRUE)
MCI_20_5 <- MutualClusteringInfo(a, f,normalize = TRUE)
MCI_20_6 <- MutualClusteringInfo(a, g,normalize = TRUE)
MCI_20_7 <- MutualClusteringInfo(a, h,normalize = TRUE)
MCI_20_8 <- MutualClusteringInfo(a, i,normalize = TRUE)
MCI_20_9 <- MutualClusteringInfo(a, j,normalize = TRUE)
MCI_20_10 <- MutualClusteringInfo(a, k,normalize = TRUE)
MCI_40_1 <- MutualClusteringInfo(a, l,normalize = TRUE)
MCI_40_2 <- MutualClusteringInfo(a, m,normalize = TRUE)
MCI_40_3 <- MutualClusteringInfo(a, n,normalize = TRUE)
MCI_40_4 <- MutualClusteringInfo(a, o,normalize = TRUE)
MCI_40_5 <- MutualClusteringInfo(a, p,normalize = TRUE)
MCI_40_6 <- MutualClusteringInfo(a, q,normalize = TRUE)
MCI_40_7 <- MutualClusteringInfo(a, r,normalize = TRUE)
MCI_40_8 <- MutualClusteringInfo(a, s,normalize = TRUE)
MCI_40_9 <- MutualClusteringInfo(a, t,normalize = TRUE)
MCI_40_10 <- MutualClusteringInfo(a, u,normalize = TRUE)
MCI_60_1 <- MutualClusteringInfo(a, v,normalize = TRUE)
MCI_60_2 <- MutualClusteringInfo(a, w,normalize = TRUE)
MCI_60_3 <- MutualClusteringInfo(a, x,normalize = TRUE)
MCI_60_4 <- MutualClusteringInfo(a, y,normalize = TRUE)
MCI_60_5 <- MutualClusteringInfo(a, z,normalize = TRUE)
MCI_60_6 <- MutualClusteringInfo(a, aa,normalize = TRUE)
MCI_60_7 <- MutualClusteringInfo(a, bb,normalize = TRUE)
MCI_60_8 <- MutualClusteringInfo(a, cc,normalize = TRUE)
MCI_60_9 <- MutualClusteringInfo(a, dd,normalize = TRUE)
MCI_60_10 <- MutualClusteringInfo(a, ee,normalize = TRUE)
MCI_80_1 <- MutualClusteringInfo(a, ff,normalize = TRUE)
MCI_80_2 <- MutualClusteringInfo(a, gg,normalize = TRUE)
MCI_80_3 <- MutualClusteringInfo(a, hh,normalize = TRUE)
MCI_80_4 <- MutualClusteringInfo(a, ii,normalize = TRUE)
MCI_80_5 <- MutualClusteringInfo(a, jj,normalize = TRUE)
MCI_80_6 <- MutualClusteringInfo(a, kk,normalize = TRUE)
MCI_80_7 <- MutualClusteringInfo(a, ll,normalize = TRUE)
MCI_80_8 <- MutualClusteringInfo(a, mm,normalize = TRUE)
MCI_80_9 <- MutualClusteringInfo(a, nn,normalize = TRUE)
MCI_80_10 <- MutualClusteringInfo(a, oo,normalize = TRUE)
#Plot MutualClusteringInfo
GMCI_1 <- c(MCI_20_1, MCI_40_1, MCI_60_1, MCI_80_1)
GMCI_2 <- c(MCI_20_2, MCI_40_2, MCI_60_2, MCI_80_2)
GMCI_3 <- c(MCI_20_3, MCI_40_3, MCI_60_3, MCI_80_3)
GMCI_4 <- c(MCI_20_4, MCI_40_4, MCI_60_4, MCI_80_4)
GMCI_5 <-  c(MCI_20_5, MCI_40_5, MCI_60_5, MCI_80_5)
GMCI_6 <- c(MCI_20_6, MCI_40_6, MCI_60_6, MCI_80_6)
GMCI_7 <- c(MCI_20_7, MCI_40_7, MCI_60_7, MCI_80_7)
GMCI_8 <- c(MCI_20_8, MCI_40_8, MCI_60_8, MCI_80_8)
GMCI_9 <- c(MCI_20_9, MCI_40_9, MCI_60_9, MCI_80_9)
GMCI_10 <- c(MCI_20_10, MCI_40_10, MCI_60_10, MCI_80_10)

df_GMCI <- data.frame(GRF,GMCI_1,GMCI_2,GMCI_3,GMCI_4,GMCI_5,GMCI_6,GMCI_7,GMCI_8,GMCI_9,GMCI_10)
ggplot(df_GMCI  , aes(GRF)) +
  xlab("The level of backbone(%)")+ylab("Mutual Clustering Infomation")+labs(title="The plot of Mutual Clustering Infomation vs. the level of backbone")+
  geom_line(aes(y=GMCI_1),
            colour="red") +
  geom_line(aes(y=GMCI_2),
            colour="green") +
  geom_line(aes(y=GMCI_3),
            colour="blue") +
  geom_line(aes(y=GMCI_4),
            colour="yellow") +
  geom_line(aes(y=GMCI_5),
            colour="orange") +
  geom_line(aes(y=GMCI_6),
            colour="black")+
  geom_line(aes(y=GMCI_7),
            colour="purple") +
  geom_line(aes(y=GMCI_8),
            colour="pink") +
  geom_line(aes(y=GMCI_9),
            colour="brown") +
  geom_line(aes(y=GMCI_10),
            colour="violet")
#Calculate NyeSimilarity
NS_20_1 <- NyeSimilarity(a, b,normalize = TRUE)
NS_20_2 <- NyeSimilarity(a, c,normalize = TRUE)
NS_20_3 <- NyeSimilarity(a, d,normalize = TRUE)
NS_20_4 <- NyeSimilarity(a, e,normalize = TRUE)
NS_20_5 <- NyeSimilarity(a, f,normalize = TRUE)
NS_20_6 <- NyeSimilarity(a, g,normalize = TRUE)
NS_20_7 <- NyeSimilarity(a, h,normalize = TRUE)
NS_20_8 <- NyeSimilarity(a, i,normalize = TRUE)
NS_20_9 <- NyeSimilarity(a, j,normalize = TRUE)
NS_20_10 <- NyeSimilarity(a, k,normalize = TRUE)
NS_40_1 <- NyeSimilarity(a, l,normalize = TRUE)
NS_40_2 <- NyeSimilarity(a, m,normalize = TRUE)
NS_40_3 <- NyeSimilarity(a, n,normalize = TRUE)
NS_40_4 <- NyeSimilarity(a, o,normalize = TRUE)
NS_40_5 <- NyeSimilarity(a, p,normalize = TRUE)
NS_40_6 <- NyeSimilarity(a, q,normalize = TRUE)
NS_40_7 <- NyeSimilarity(a, r,normalize = TRUE)
NS_40_8 <- NyeSimilarity(a, s,normalize = TRUE)
NS_40_9 <- NyeSimilarity(a, t,normalize = TRUE)
NS_40_10 <- NyeSimilarity(a, u,normalize = TRUE)
NS_60_1 <- NyeSimilarity(a, v,normalize = TRUE)
NS_60_2 <- NyeSimilarity(a, w,normalize = TRUE)
NS_60_3 <- NyeSimilarity(a, x,normalize = TRUE)
NS_60_4 <- NyeSimilarity(a, y,normalize = TRUE)
NS_60_5 <- NyeSimilarity(a, z,normalize = TRUE)
NS_60_6 <- NyeSimilarity(a, aa,normalize = TRUE)
NS_60_7 <- NyeSimilarity(a, bb,normalize = TRUE)
NS_60_8 <- NyeSimilarity(a, cc,normalize = TRUE)
NS_60_9 <- NyeSimilarity(a, dd,normalize = TRUE)
NS_60_10 <- NyeSimilarity(a, ee,normalize = TRUE)
NS_80_1 <- NyeSimilarity(a, ff,normalize = TRUE)
NS_80_2 <- NyeSimilarity(a, gg,normalize = TRUE)
NS_80_3 <- NyeSimilarity(a, hh,normalize = TRUE)
NS_80_4 <- NyeSimilarity(a, ii,normalize = TRUE)
NS_80_5 <- NyeSimilarity(a, jj,normalize = TRUE)
NS_80_6 <- NyeSimilarity(a, kk,normalize = TRUE)
NS_80_7 <- NyeSimilarity(a, ll,normalize = TRUE)
NS_80_8 <- NyeSimilarity(a, mm,normalize = TRUE)
NS_80_9 <- NyeSimilarity(a, nn,normalize = TRUE)
NS_80_10 <- NyeSimilarity(a, oo,normalize = TRUE)
#Plot NyeSimilarity
GNS_1 <- c(NS_20_1, NS_40_1, NS_60_1, NS_80_1)
GNS_2 <- c(NS_20_2, NS_40_2, NS_60_2, NS_80_2)
GNS_3 <- c(NS_20_3, NS_40_3, NS_60_3, NS_80_3)
GNS_4 <- c(NS_20_4, NS_40_4, NS_60_4, NS_80_4)
GNS_5 <-  c(NS_20_5, NS_40_5, NS_60_5, NS_80_5)
GNS_6 <- c(NS_20_6, NS_40_6, NS_60_6, NS_80_6)
GNS_7 <- c(NS_20_7, NS_40_7, NS_60_7, NS_80_7)
GNS_8 <- c(NS_20_8, NS_40_8, NS_60_8, NS_80_8)
GNS_9 <- c(NS_20_9, NS_40_9, NS_60_9, NS_80_9)
GNS_10 <- c(NS_20_10, NS_40_10, NS_60_10, NS_80_10)

df_GNS <- data.frame(GRF,GNS_1,GNS_2,GNS_3,GNS_4,GNS_5,GNS_6,GNS_7,GNS_8,GNS_9,GNS_10)
ggplot(df_GNS  , aes(GRF)) +
  xlab("The level of backbone(%)")+ylab("Nye Similarity")+labs(title="The plot of Nye Similarity vs. the level of backbone")+
  geom_line(aes(y=GNS_1),
            colour="red") +
  geom_line(aes(y=GNS_2),
            colour="green") +
  geom_line(aes(y=GNS_3),
            colour="blue") +
  geom_line(aes(y=GNS_4),
            colour="yellow") +
  geom_line(aes(y=GNS_5),
            colour="orange") +
  geom_line(aes(y=GNS_6),
            colour="black")+
  geom_line(aes(y=GNS_7),
            colour="purple") +
  geom_line(aes(y=GNS_8),
            colour="pink") +
  geom_line(aes(y=GNS_9),
            colour="brown") +
  geom_line(aes(y=GNS_10),
            colour="violet")
#Calculate JaccardRobinsonFoulds
JRF_20_1 <- JaccardRobinsonFoulds(a, b,normalize = TRUE)
JRF_20_2 <- JaccardRobinsonFoulds(a, c,normalize = TRUE)
JRF_20_3 <- JaccardRobinsonFoulds(a, d,normalize = TRUE)
JRF_20_4 <- JaccardRobinsonFoulds(a, e,normalize = TRUE)
JRF_20_5 <- JaccardRobinsonFoulds(a, f,normalize = TRUE)
JRF_20_6 <- JaccardRobinsonFoulds(a, g,normalize = TRUE)
JRF_20_7 <- JaccardRobinsonFoulds(a, h,normalize = TRUE)
JRF_20_8 <- JaccardRobinsonFoulds(a, i,normalize = TRUE)
JRF_20_9 <- JaccardRobinsonFoulds(a, j,normalize = TRUE)
JRF_20_10 <- JaccardRobinsonFoulds(a, k,normalize = TRUE)
JRF_40_1 <- JaccardRobinsonFoulds(a, l,normalize = TRUE)
JRF_40_2 <- JaccardRobinsonFoulds(a, m,normalize = TRUE)
JRF_40_3 <- JaccardRobinsonFoulds(a, n,normalize = TRUE)
JRF_40_4 <- JaccardRobinsonFoulds(a, o,normalize = TRUE)
JRF_40_5 <- JaccardRobinsonFoulds(a, p,normalize = TRUE)
JRF_40_6 <- JaccardRobinsonFoulds(a, q,normalize = TRUE)
JRF_40_7 <- JaccardRobinsonFoulds(a, r,normalize = TRUE)
JRF_40_8 <- JaccardRobinsonFoulds(a, s,normalize = TRUE)
JRF_40_9 <- JaccardRobinsonFoulds(a, t,normalize = TRUE)
JRF_40_10 <- JaccardRobinsonFoulds(a, u,normalize = TRUE)
JRF_60_1 <- JaccardRobinsonFoulds(a, v,normalize = TRUE)
JRF_60_2 <- JaccardRobinsonFoulds(a, w,normalize = TRUE)
JRF_60_3 <- JaccardRobinsonFoulds(a, x,normalize = TRUE)
JRF_60_4 <- JaccardRobinsonFoulds(a, y,normalize = TRUE)
JRF_60_5 <- JaccardRobinsonFoulds(a, z,normalize = TRUE)
JRF_60_6 <- JaccardRobinsonFoulds(a, aa,normalize = TRUE)
JRF_60_7 <- JaccardRobinsonFoulds(a, bb,normalize = TRUE)
JRF_60_8 <- JaccardRobinsonFoulds(a, cc,normalize = TRUE)
JRF_60_9 <- JaccardRobinsonFoulds(a, dd,normalize = TRUE)
JRF_60_10 <- JaccardRobinsonFoulds(a, ee,normalize = TRUE)
JRF_80_1 <- JaccardRobinsonFoulds(a, ff,normalize = TRUE)
JRF_80_2 <- JaccardRobinsonFoulds(a, gg,normalize = TRUE)
JRF_80_3 <- JaccardRobinsonFoulds(a, hh,normalize = TRUE)
JRF_80_4 <- JaccardRobinsonFoulds(a, ii,normalize = TRUE)
JRF_80_5 <- JaccardRobinsonFoulds(a, jj,normalize = TRUE)
JRF_80_6 <- JaccardRobinsonFoulds(a, kk,normalize = TRUE)
JRF_80_7 <- JaccardRobinsonFoulds(a, ll,normalize = TRUE)
JRF_80_8 <- JaccardRobinsonFoulds(a, mm,normalize = TRUE)
JRF_80_9 <- JaccardRobinsonFoulds(a, nn,normalize = TRUE)
JRF_80_10 <- JaccardRobinsonFoulds(a, oo,normalize = TRUE)
#Plot JaccardRobinsonFoulds
GJRF_1 <- c(JRF_20_1, JRF_40_1, JRF_60_1, JRF_80_1)
GJRF_2 <- c(JRF_20_2, JRF_40_2, JRF_60_2, JRF_80_2)
GJRF_3 <- c(JRF_20_3, JRF_40_3, JRF_60_3, JRF_80_3)
GJRF_4 <- c(JRF_20_4, JRF_40_4, JRF_60_4, JRF_80_4)
GJRF_5 <-  c(JRF_20_5, JRF_40_5, JRF_60_5, JRF_80_5)
GJRF_6 <- c(JRF_20_6, JRF_40_6, JRF_60_6, JRF_80_6)
GJRF_7 <- c(JRF_20_7, JRF_40_7, JRF_60_7, JRF_80_7)
GJRF_8 <- c(JRF_20_8, JRF_40_8, JRF_60_8, JRF_80_8)
GJRF_9 <- c(JRF_20_9, JRF_40_9, JRF_60_9, JRF_80_9)
GJRF_10 <- c(JRF_20_10, JRF_40_10, JRF_60_10, JRF_80_10)

df_GJRF <- data.frame(GRF,GJRF_1,GJRF_2,GJRF_3,GJRF_4,GJRF_5,GJRF_6,GJRF_7,GJRF_8,GJRF_9,GJRF_10)
ggplot(df_GJRF  , aes(GRF)) + scale_y_reverse() +
  xlab("The level of backbone(%)")+ylab("Jaccard Robinson Foulds")+labs(title="The plot of Jaccard Robinson Foulds vs. the level of backbone")+
  geom_line(aes(y=GJRF_1),
            colour="red") +
  geom_line(aes(y=GJRF_2),
            colour="green") +
  geom_line(aes(y=GJRF_3),
            colour="blue") +
  geom_line(aes(y=GJRF_4),
            colour="yellow") +
  geom_line(aes(y=GJRF_5),
            colour="orange") +
  geom_line(aes(y=GJRF_6),
            colour="black")+
  geom_line(aes(y=GJRF_7),
            colour="purple") +
  geom_line(aes(y=GJRF_8),
            colour="pink") +
  geom_line(aes(y=GJRF_9),
            colour="brown") +
  geom_line(aes(y=GJRF_10),
            colour="violet")
#Calculate MatchingSplitInfoDistance
MSD_20_1 <- MatchingSplitInfoDistance(a, b,normalize = TRUE)
MSD_20_2 <- MatchingSplitInfoDistance(a, c,normalize = TRUE)
MSD_20_3 <- MatchingSplitInfoDistance(a, d,normalize = TRUE)
MSD_20_4 <- MatchingSplitInfoDistance(a, e,normalize = TRUE)
MSD_20_5 <- MatchingSplitInfoDistance(a, f,normalize = TRUE)
MSD_20_6 <- MatchingSplitInfoDistance(a, g,normalize = TRUE)
MSD_20_7 <- MatchingSplitInfoDistance(a, h,normalize = TRUE)
MSD_20_8 <- MatchingSplitInfoDistance(a, i,normalize = TRUE)
MSD_20_9 <- MatchingSplitInfoDistance(a, j,normalize = TRUE)
MSD_20_10 <- MatchingSplitInfoDistance(a, k,normalize = TRUE)
MSD_40_1 <- MatchingSplitInfoDistance(a, l,normalize = TRUE)
MSD_40_2 <- MatchingSplitInfoDistance(a, m,normalize = TRUE)
MSD_40_3 <- MatchingSplitInfoDistance(a, n,normalize = TRUE)
MSD_40_4 <- MatchingSplitInfoDistance(a, o,normalize = TRUE)
MSD_40_5 <- MatchingSplitInfoDistance(a, p,normalize = TRUE)
MSD_40_6 <- MatchingSplitInfoDistance(a, q,normalize = TRUE)
MSD_40_7 <- MatchingSplitInfoDistance(a, r,normalize = TRUE)
MSD_40_8 <- MatchingSplitInfoDistance(a, s,normalize = TRUE)
MSD_40_9 <- MatchingSplitInfoDistance(a, t,normalize = TRUE)
MSD_40_10 <- MatchingSplitInfoDistance(a, u,normalize = TRUE)
MSD_60_1 <- MatchingSplitInfoDistance(a, v,normalize = TRUE)
MSD_60_2 <- MatchingSplitInfoDistance(a, w,normalize = TRUE)
MSD_60_3 <- MatchingSplitInfoDistance(a, x,normalize = TRUE)
MSD_60_4 <- MatchingSplitInfoDistance(a, y,normalize = TRUE)
MSD_60_5 <- MatchingSplitInfoDistance(a, z,normalize = TRUE)
MSD_60_6 <- MatchingSplitInfoDistance(a, aa,normalize = TRUE)
MSD_60_7 <- MatchingSplitInfoDistance(a, bb,normalize = TRUE)
MSD_60_8 <- MatchingSplitInfoDistance(a, cc,normalize = TRUE)
MSD_60_9 <- MatchingSplitInfoDistance(a, dd,normalize = TRUE)
MSD_60_10 <- MatchingSplitInfoDistance(a, ee,normalize = TRUE)
MSD_80_1 <- MatchingSplitInfoDistance(a, ff,normalize = TRUE)
MSD_80_2 <- MatchingSplitInfoDistance(a, gg,normalize = TRUE)
MSD_80_3 <- MatchingSplitInfoDistance(a, hh,normalize = TRUE)
MSD_80_4 <- MatchingSplitInfoDistance(a, ii,normalize = TRUE)
MSD_80_5 <- MatchingSplitInfoDistance(a, jj,normalize = TRUE)
MSD_80_6 <- MatchingSplitInfoDistance(a, kk,normalize = TRUE)
MSD_80_7 <- MatchingSplitInfoDistance(a, ll,normalize = TRUE)
MSD_80_8 <- MatchingSplitInfoDistance(a, mm,normalize = TRUE)
MSD_80_9 <- MatchingSplitInfoDistance(a, nn,normalize = TRUE)
MSD_80_10 <- MatchingSplitInfoDistance(a, oo,normalize = TRUE)
#Plot MatchingSplitInfoDistance
GMSD_1 <- c(MSD_20_1, MSD_40_1, MSD_60_1, MSD_80_1)
GMSD_2 <- c(MSD_20_2, MSD_40_2, MSD_60_2, MSD_80_2)
GMSD_3 <- c(MSD_20_3, MSD_40_3, MSD_60_3, MSD_80_3)
GMSD_4 <- c(MSD_20_4, MSD_40_4, MSD_60_4, MSD_80_4)
GMSD_5 <- c(MSD_20_5, MSD_40_5, MSD_60_5, MSD_80_5)
GMSD_6 <- c(MSD_20_6, MSD_40_6, MSD_60_6, MSD_80_6)
GMSD_7 <- c(MSD_20_7, MSD_40_7, MSD_60_7, MSD_80_7)
GMSD_8 <- c(MSD_20_8, MSD_40_8, MSD_60_8, MSD_80_8)
GMSD_9 <- c(MSD_20_9, MSD_40_9, MSD_60_9, MSD_80_9)
GMSD_10 <- c(MSD_20_10, MSD_40_10, MSD_60_10, MSD_80_10)

df_GMSD <- data.frame(GRF,GMSD_1,GMSD_2,GMSD_3,GMSD_4,GMSD_5,GMSD_6,GMSD_7,GMSD_8,GMSD_9,GMSD_10)
ggplot(df_GMSD  , aes(GRF)) + scale_y_reverse() +
  xlab("The level of backbone(%)")+ylab("Matching Split Info Distance")+labs(title="The plot of Matching Split Info Distance vs. the level of backbone")+
  geom_line(aes(y=GMSD_1),
            colour="red") +
  geom_line(aes(y=GMSD_2),
            colour="green") +
  geom_line(aes(y=GMSD_3),
            colour="blue") +
  geom_line(aes(y=GMSD_4),
            colour="yellow") +
  geom_line(aes(y=GMSD_5),
            colour="orange") +
  geom_line(aes(y=GMSD_6),
            colour="black")+
  geom_line(aes(y=GMSD_7),
            colour="purple") +
  geom_line(aes(y=GMSD_8),
            colour="pink") +
  geom_line(aes(y=GMSD_9),
            colour="brown") +
  geom_line(aes(y=GMSD_10),
            colour="violet")




