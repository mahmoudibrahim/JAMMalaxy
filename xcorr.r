###############################
#### Fragment Length Calculator
#### R script		
###############################


# ========================= 
# Required Libraries check 
# ========================= 
if ((is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package'parallel' is required. Please install it!")
}
suppressPackageStartupMessages(library("parallel"))
#=======================> DONE! 




# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
ibed = NA # input bed file
sFile = NA # chromosome size
storeFile = NA # file to store result

#Parsing arguments and storing values
for (each.arg in args) {
	#input bed file
	if (grepl('^-ibed=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				ibed <- arg.split[2]
		} else {
			stop('No input bed file')
		} 
	}
	if (grepl('^-s=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			sFile <- arg.split[2] 
		} else {
			stop('No chromosome size file')
		}
	}
	if (grepl('^-rl=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			rl <- arg.split[2] 
		} else {
			stop('Read length missing')
		}
	}
	if (grepl('^-d=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			storeFile <- arg.split[2] 
		} else {
			stop('No file to store result')
		}
	}
	if (grepl('^-p=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			cornum <- as.numeric(arg.split[2]) 
		} else {
			stop('No number of cores given')
		}
	}
}

#Read in variables
chromosomes = read.table(sFile, header=FALSE)
chromName = chromosomes$V1 #which chromosome
chromSize = as.numeric(chromosomes$V2) #chromosome size
rm(chromosomes)

ibed = strsplit(ibed, ",", fixed = TRUE)[[1]]
rl = as.numeric(strsplit(rl, ",", fixed = TRUE)[[1]])
#=======================> DONE! 


counting = function(readst, br){
	return(hist(readst, breaks = br, plot = FALSE)$counts)
}


countingreads = function(bedfile, readlen, filelist) {
	reads = read.table(bedfile, header = FALSE)
	reads = cbind(as.character(reads[[1]]), as.character(reads[[2]]))
	o = which(filelist == bedfile)
	tri = list(as.numeric(reads[reads[,2] == "+",,drop = FALSE][,1]), (as.numeric(reads[reads[,2] == "-",,drop = FALSE][,1])) + readlen[o] - 1)
	readcounts = mclapply(tri, counting, br = genomevec, mc.cores = cornum)
	return(readcounts)
}

xc = function(countlist) {
	crossCorrelation = ccf(countlist[[2]], countlist[[1]], plot = FALSE); #xcorr
	crossCorrelation$lag = crossCorrelation$lag * 20; #correct lag for counts window
	maxCorr = which.max(crossCorrelation$acf);
	maxCorr = abs(crossCorrelation$lag[maxCorr]);
	return(maxCorr)
}


# ================== 
# Read Start Count
# ==================
#make the window vector (a vector with start positions of chromosome windows)
genomevec = seq(0, chromSize, by = 20);
if (max(genomevec) < chromSize) {
	genomevec = append(genomevec, chromSize);
}
datain = mclapply(ibed, countingreads, readlen = rl, filelist = ibed, mc.cores = cornum)
#=======================> DONE! 


# ================== 
# Cross Correlation
# ==================
xcorr = unlist(mclapply(datain, xc, mc.cores = cornum))
#=======================> DONE! 


# ===================== 
# Write result to File
# =====================
for (i in 1:length(xcorr)) {
	filename = strsplit(ibed[i], "/", fixed = TRUE)[[1]]
	filename = filename[length(filename)]
	filename = strsplit(filename, ".", fixed = TRUE)[[1]]
	filename = filename[3]
	if ((xcorr[i] <= 500) && (xcorr[i] >= 50)) { #only write this to file if xcorr value was plausible
		message(paste0(chromName, ", ", filename, ": Ok!"))
		write(paste(chromName, xcorr[i], sep = "\t"), file = paste0(storeFile, "/xc.", filename, ".tab"), append = TRUE)
	} else {
		message(paste0(chromName, ", ", filename, ": Value Not Used!"))
	}
}
#=======================> DONE!
