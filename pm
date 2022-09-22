#!/usr/bin/env Rscript
bindir <- ""
cseries.dir <- ""
functions.dir <- ""
##------------------------------------------------------------------------------

require("ape",quietly=TRUE)
require("phangorn",quietly=TRUE)
require("quadprog",quietly=TRUE)

# Make sure mammal_functions and pm_functions are in the same directory as this program
source(paste0(functions.dir,"mammal_functions")) 
source(paste0(functions.dir,"pm_functions"))

##------------------------------------------------------------------------------
## PM (MAMMaL Extension for Entropy and Rates) - adds Entropy, Upper/Lower partitions, Invariable flag, and includes low rates

## External programs
dgpe <- paste(bindir,"dgpe",sep="")
mult.data <- paste(bindir,"mult-data",sep="")

FnameE <- function(fname){
  if(!file.exists(fname)) stop(paste(fname,"does not exist"))
  return(fname)
}

## Process Command Line
args <- commandArgs()
iarg <- length(args)

nclass <- c(NULL,NULL) # number of classes ([1] = high partition, [2] = low partition)
iqtreefile <- seqfile <- treefile <- NULL
outfile <- "esmodel" # name of output .nex file
q <- 0.75 # quantile for rate estimation
lwt.needed <- output <- TRUE
start.frs.type <- "hclust"
log.method <- "both"
invar <- add.invar.class <- plusF <- cluster.set <- sort <- epsilon <- FALSE
saveFrequencyFiles <- savePartitionFiles <- saveTempFiles <- useSuffix <- FALSE
partition.mode <- "E" # partition by entropy by default
C <- 1.0e-10 # penalty parameter
min.lwt <- 1.0e-10 # minimum likelihood weight

while(iarg>=6){
  not.an.option <- TRUE
  opt <- args[iarg]
  val <- NULL
  if(substring(opt,1,1)!='-'){
    if(iarg == 6) stop(paste(opt,"is not an option"))
    val <- opt
    iarg <- iarg - 1
    opt <- args[iarg]
  }
  if(opt=="-s"){
    if(is.null(val)) stop("seqfile not specified in -s seqfile")
    seqfile <- orig.seqfile <- FnameE(val); not.an.option <- FALSE
  }
  if(opt=="-t"){
    if(is.null(val)) stop("treefile not specified in -t treefile")
    treefile <- FnameE(val); not.an.option <- FALSE
  }
  if(opt=="-o"){
    if(is.null(val)) stop("output file name not specified in -o outfile")
    outfile <- val; not.an.option <- FALSE
  }
  if(opt=="-c"){
    if(is.null(val)) stop("number of classes not specified in -c nclasses")
    nclass[1] <- nclass[2] <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-ch"){
    if(is.null(val)) stop("number of classes in high partition not specified in -ch nhigh")
    nclass[1] <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-cl"){
    if(is.null(val)) stop("number of classes in low partition not specified in -cl nlow")
    nclass[2] <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-m"){
    output <- FALSE; not.an.option <- FALSE
  }
  if(opt=="-l"){
    lwt.needed <- FALSE; not.an.option <- FALSE
  }
  if(opt=="-lm"){
    min.lwt <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-q"){
    if(is.null(val)) stop("quantile not specified in -q quantile")
    q <- as.numeric(val); not.an.option <- FALSE
    if(q == 0) nclass[2] <- 0
  }
  if(opt=="-C"){
    if(is.null(val)) stop("penalty not specified in -C penalty")
    C <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-d"){
    plusF <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-p"){
    if(is.null(val)) stop("partition mode not specified in -p mode")
    val <- toupper(val)
    if(!val %in% c("E","R","K")) stop(paste0("invalid partition mode (",val,"); valid modes are: E (entropy), R (rates), K (k_eff)"))
    partition.mode <- val; not.an.option <- FALSE
    
  }
  if(opt=="-f"){
    if(is.null(val)) stop("cluster type not specified in -f type (accepted values: C (c-series), H (hclust) )")
    val <- toupper(val)
    if(!val %in% c("C","H")) stop(paste0("invalid cluster type (",val,"); valid types are: C (c-series), H (hclust)"))
    if(val == "C") start.frs.type <- "c-series"
    if(val == "H") start.frs.type <- "hclust"
    not.an.option <- FALSE
    cluster.set <- TRUE
  }
  if(opt=="-log"){
    if(is.null(val)) stop("logging method not specified in -log method (accepted values: N (no log), F (log file), X (log in nexus file) )")
    val <- toupper(val)
    if(!val %in% c("N","F","X","B")) stop(paste0("invalid logging method (",val,"); valid types are: N (no log), F (log file), X (log in nexus file)"))
    if(val == "N") log.method <- "none"
    if(val == "F") log.method <- "file"
    if(val == "X") log.method <- "nexus"
    if(val == "B") log.method <- "both"
    not.an.option <- FALSE
  }
  if(opt=="-I"){
    add.invar.class <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-ri"){
    invar <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-pf"){
    savePartitionFiles <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-ff"){
    saveFrequencyFiles <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-tmp"){
    saveTempFiles <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-suf"){
    useSuffix <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-sort"){
    sort <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-e"){
    epsilon <- TRUE; not.an.option <- FALSE
  }
  if(not.an.option) stop(paste(opt,"is not an option"))
  iarg <- iarg - 1 
}

if(is.null(seqfile)) stop("sequence file needed: -s seqfile")
if(is.null(nclass[1]) | is.null(nclass[2])) stop("must specify number of classes in upper and lower partitions: -cl nlow; -ch nhigh (alternatively: -c nclasses [for both])")

# start timer
T1 <- Sys.time()

# determine if we should use c-series or not
if(!cluster.set & partition.mode == "R" & nclass[1] %in% seq(0,60,10) & nclass[2] %in% seq(0,60,10)) start.frs.type <- "c-series"

# Keep track of events to log at the end
file.log <- c(as.character(Sys.time()))
Log <- function(line) {
  cat(paste0(line,"\n"))
  if(grepl("\n",line)) line <- unlist(strsplit(line,"\n"))
  file.log <<- append(file.log,line)
}

# Generate a random file suffix to support multiple threads at once
if(useSuffix) {
  file.suffix <- as.character(round(runif(1) * 100000)) 
  while(nchar(file.suffix) < 4) file.suffix <- as.character(round(runif(1) * 100000)) # use at least a 4-digit number
  Log(paste("Generated file suffix:",file.suffix))
  file.suffix <- paste0("_",file.suffix) #use underscore in front for more readable file names
} else file.suffix <- ""
# Utility function to imitate paste0 and add the file suffix to the end
file.format <- function(...) { 
  fname <- paste0(...,collapse="")
  if(grepl(file.suffix,fname)){
    # remove the suffix if it already exists so we can append it back to the end
    # assumption: the suffix only appears once in the file name
    fname <- removeSuffix(fname,file.suffix) 
  }
  return(paste0(fname,file.suffix)) 
} 

# Check if the treefile has any negatives and zero them
if(partition.mode == "R") {
  if(is.null(treefile)) {
    generateTree(seqfile) # generate a tree file if using rates and one not provided
    treefile <- paste0(seqfile,".tree")
    Log("No tree file provided - generating new neighbour-joined tree using the LG model")
  }
  valid <- verifyTree(treefile,suffix=file.suffix)
  if(!valid){
    treefile <- file.format(treefile,".tmp") # autogenerated by verifyTree if false
    Log("WARNING: Tree file contains invalid edge lengths. Proceeding with these edges zeroed.")
  }
}

total.classes <- (nclass[1]+ifelse(nclass[1]>0,ifelse(plusF,1,0),0))+(nclass[2]+ifelse(nclass[2]>0,ifelse(plusF,1,0),0))+ifelse(add.invar.class,1,0)
Log(paste("Running PM with the following parameters:\n",
            "Sequence file: ",seqfile,"\t",ifelse(is.null(treefile),"",paste("Tree file: ",treefile)),"\n",
            "Output file: ",paste0(outfile,".nex (Total number of classes: ",total.classes,")"),"\n",
            "Partition mode: ",ifelse(partition.mode=="R","Rates",""),ifelse(partition.mode=="E","Entropy",""),ifelse(partition.mode=="K","K_eff",""),"\n",
            "High partition classes: ",nclass[1],"\tLow partition classes: ",nclass[2],"\n",
            "Quantile: ",q,"\tCluster type: ",start.frs.type,"\n",
            ifelse(invar,"Remove invariant: TRUE\n",""),
            ifelse(add.invar.class,"Add invariant class: TRUE\n",""),
            ifelse(sort,"Sort by Entropy: TRUE\n",""),
            ifelse(epsilon,"Add Epsilon to H-Clust: TRUE","")
            )
      )

# Count and output number of sequences and sites in the sequence file
count <- seqCount(seqfile) 
Log(paste("Number of sequences:",count[1],";","Number of sites:",count[2]))

cat("Processing sequence file...\n")

# Remove any empty sites from the sequence
num_e <- removeEmptySites(seqfile,file.suffix)
if(num_e[1] > 0) {
  Log(paste("Removed",num_e[1],"empty sites (new count:",num_e[2],")"))
  seqfile <- file.format(seqfile,".e")
}

# Detect invariant sites and use new file with them removed if the flag is set
invar.results <- detectInvariant(seqfile,invar,suffix=file.suffix)
num_i <- c(invar.results[[1]],invar.results[[2]])
invar.class <- invar.results[[3]]
if(invar) {
  Log(paste("Removed",num_i[1],"invariant sites (new count:",num_i[2],")"))
  if(num_e[1] > 0) did.remove <- file.remove(seqfile) # remove temp sequence file if one was generated
  seqfile <- file.format(seqfile,".i")
} else Log(paste("Detected",num_i[1],"invariant sites"))

# Try to sort the sequence file so we get more consistent output from different permutations
if(sort) {
  sortSequenceFile(seqfile,suffix=file.suffix)
  sorted.seq <- file.format(seqfile,".s")
  if(num_e[1] > 0 | invar) did.remove <- file.remove(seqfile) # remove temp sequence file if one was generated
  seqfile <- sorted.seq
  Log(paste("Sorted sequence file by Site Entropy"))
}

# We should get two files, one for high partition and one for low
expectedFiles <- vector(length=2) 

## Extract high and low rate/entropy/k_eff sites
cat("Creating high and low partitions...\n")
if(partition.mode == "R") {
  # Run dgpe to get the high rate sites
  system(paste(dgpe,"-i",seqfile,"-t",treefile,"-o",file.format(seqfile,".highrate"),"-q",q," >",file.format("tmp.out")))
  # Use output from dgpe to get low rate sites
  extractLowRates(seqfile,"rate_est.dat",q,suffix=file.suffix) 
  expectedFiles[1] <- file.format(seqfile,".highrate")
  expectedFiles[2] <- file.format(seqfile,".lowrate")
} else if(partition.mode == "E") {
  partition(seqfile,func = getEntropy, suffix=file.suffix)
  expectedFiles[1] <- file.format(seqfile,".highentropy")
  expectedFiles[2] <- file.format(seqfile,".lowentropy")
} else if(partition.mode == "K") {
  partition(seqfile,func = getEffectiveK, suffix=file.suffix)
  expectedFiles[1] <- file.format(seqfile,".highkeff")
  expectedFiles[2] <- file.format(seqfile,".lowkeff")
}

ntaxa <- vector(length=2)
lwt <- list(NULL,NULL)

## Likelihood weights
cat("Calculating likelihood weights...\n")
for(i in 1:2) {
  
  if(nclass[i] == 0) {next} # skip if there are no classes for this partition
  
  ntaxa[i] <- as.numeric(scan(expectedFiles[i],n=1,quiet=TRUE))
  if(lwt.needed){
    lwt[[i]] <- CalculateWeightsE(seqfile=expectedFiles[i],bindir=bindir,clean=FALSE,suffix=file.suffix)$w
  }else{
    lwt[[i]] <- rep(1,ntaxa[i])
  }
  
  ## ## Use this to allow 0 weights
  ## if(sum(lwt<0)>0) lwt[lwt<0] <- 0
  ## lwt <- ntaxa*lwt/sum(lwt)
  ## Minimum weight >= min.lwt
  zeroWeights <- function(lw,nt)
  {
    if(sum(lw<=min.lwt)>0){
      lw[lw<=min.lwt] <- 0
      a <- nt*(1-min.lwt)/sum(lw)
      return(a*lw + min.lwt)
    }
  }
  
  lwt[[i]] <- zeroWeights(lwt[[i]],ntaxa[i])
  write(format(lwt[[i]],sci=TRUE,digits=16),ncol=1,file=file.format("tmp.lwt",i))
}

frs <- list(NULL,NULL)

# Do the following for both high (i=1) and low (i=2) partitions
cat("Processing data...\n")
for(i in 1:2) {
  
  if(nclass[i] == 0) {next} # skip if there are no classes for this partition
  
  ## Starting frequencies
  if(start.frs.type =="hclust"){
    frs[[i]] <- HclustCenters(DataFrequencies(expectedFiles[i],
                                         clean=TRUE,bindir=bindir,suffix=file.suffix),
                         hclust.type="hclust",
                         nclass=nclass[i],dmethod="manhattan")
    if(epsilon){ frs[[i]] <- addEpsilonToHclust(frs[[i]]) }
  }
  if(start.frs.type =="c-series"){
    frs[[i]] <- scan(paste(cseries.dir,"C",nclass[i],".aafreq.dat",sep=""),quiet=TRUE)
    frs[[i]] <- matrix(frs[[i]],ncol=20,byrow=TRUE)
  }
  write(format(t(frs[[i]]),sci=TRUE,digits=16),file=file.format("tmp.frs",i),ncol=20)
  if(plusF){
    system(paste0(bindir,"charfreq 20 < ",seqfile," >> tmp.frs",i,file.suffix))
    nclass[i] <- nclass[i]+1
  }

  ## mult-mix-lwt
  cmdline <- paste(bindir,"mult-mix-lwt -i ",expectedFiles[i],
                   paste0(" -l tmp.lwt",i,file.suffix),
                   # ifelse(invar," -0",""),
                   #" -c ",ifelse(invar,nclass[i]+1,nclass[i]),
                   " -c ",nclass[i],
                   paste0(" -f tmp.frs",i,file.suffix),
                   paste0(" -o estimated-frequencies",i,file.suffix),
                   " -w estimated-weights",
                   ifelse(plusF," -d ",""),
                   " -p ",
                   ifelse(C>0,paste(" -C",C),""),                 
                   " > tmp.err",file.suffix,sep="")
  system(cmdline)
}

lnl <-  scan(file.format("tmp.err"),quiet=TRUE)

# stop timer
T2 <- Sys.time()
seconds <- round(as.numeric(difftime(T2, T1, units = "secs")))
minutes <- hours <- 0
if(seconds > 60) { minutes <- floor(seconds / 60); seconds <- seconds - (minutes * 60) }
if(minutes > 60) { hours <- floor(minutes / 60); minutes <- minutes - (hours * 60) }
timestring <- paste(sprintf("%02d",hours),sprintf("%02d",minutes),sprintf("%02d",seconds),sep=":")
file.log[1] <- paste0(file.log[1],"     Execution time: ",timestring)

## Output
if(output){
  cat("Generating output...\n")
  highExists <- (nclass[1] > 0)
  lowExists <- (nclass[2] > 0)
  if(highExists) {
    fr1 <- matrix(scan(file.format("estimated-frequencies1"),quiet=TRUE),ncol=20,byrow=TRUE)
    CreateIQFreqfile(fr1,file.format(outfile,".high.nex"))
    if(!lowExists){
      if(add.invar.class) addInvariantClassToNexus(paste0(outfile,".high.nex"), invar.class)
      if(log.method == "nexus" | log.method == "both") nexusHeader(file.format(outfile,".high.nex"), file.log)
    }
  }
  if(lowExists) {
    fr2 <- matrix(scan(file.format("estimated-frequencies2"),quiet=TRUE),ncol=20,byrow=TRUE)
    CreateIQFreqfile(fr2,file.format(outfile,".low.nex"))
    if(!highExists){
      if(add.invar.class) addInvariantClassToNexus(paste0(outfile,".low.nex"), invar.class)
      if(log.method == "nexus" | log.method == "both") nexusHeader(file.format(outfile,".low.nex"), file.log)
    }
  }
  
  # merge nexus files and delete them if they exist
  if(highExists & lowExists) {
    mergeNexFiles(file.format(outfile,".high.nex"),file.format(outfile,".low.nex"),outfile)
    did.remove <- file.remove(c(file.format(outfile,".high.nex"),file.format(outfile,".low.nex")))
    if(add.invar.class) addInvariantClassToNexus(paste0(outfile,".nex"), invar.class)
    if(log.method == "nexus" | log.method == "both") nexusHeader(paste0(outfile,".nex"), file.log)
  }
  # if they don't both exist, rename the .high or .low file to the appropriate output
  else {
    if(highExists) rn <- file.rename(file.format(outfile,".high.nex"), file.format(outfile,".nex"))
    if(lowExists)  rn <- file.rename(file.format(outfile,".low.nex"), file.format(outfile,".nex"))
  }
}

# Write log file
if(log.method == "file" | log.method == "both") {
  fname <- paste0("pm_output",file.suffix,".log")
  if(file.exists(fname)){
    logfile <- file(fname)
    oldLines <- readLines(fname)
    writeLines(c(file.log,"-------------------------------",oldLines),logfile)
  }else{
    logfile <- file(fname)
    writeLines(file.log,logfile)
  }
  close(logfile)
}

## Remove all temporary files
tmp.files <- c("tmp.out", "tmp.Sigma",
               "tmp.frs1","tmp.frs2","tmp.err",
               "tmp.lwt1","tmp.lwt2")
tmp.files <- unlist(lapply(tmp.files,file.format))
tmp.files <- append(tmp.files, c("rate_est.dat","estimated-weights"))
if(partition.mode == "R") { if(grepl(".tmp",treefile)) tmp.files <- append(tmp.files,treefile) }
if(saveTempFiles) tmp.files <- vector() # keep the temp files by resetting tmp.files to empty vector
if(seqfile != orig.seqfile) tmp.files <- append(tmp.files,seqfile)
if(!saveFrequencyFiles) tmp.files <- append(tmp.files, 
             c(file.format("estimated-frequencies1"), file.format("estimated-frequencies2")) )
if(!savePartitionFiles) tmp.files <- append(tmp.files,expectedFiles)

to.remove <- vector()
for(fname in tmp.files) {
  if(file.exists(fname)) to.remove <- append(to.remove,fname)
}
did.remove <- file.remove(to.remove)

cat("Done\n")
