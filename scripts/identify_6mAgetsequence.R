calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  if(drawPlot){  #if TRUE, draw the plot
    plot(1:length(inputVector), inputVector,type="l",...)
    b <- y_cutoff-(slope* xPt)
    abline(v= xPt,h= y_cutoff,lty=2,col=8)
    points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
    abline(coef=c(b,slope),col=2)
    title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
  }  
  return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}
convert_stitched_to_bed <- function(inputStitched, trackName, trackDescription, outputFile, 
                                    splitSuper = TRUE, score = c(), superRows = c(), 
                                    baseColor = "0,0,0", superColor = "255,0,0") {
  outMatrix <- matrix(data = "", ncol = 5, nrow = nrow(inputStitched))  
  
  outMatrix[, 1] <- as.character(inputStitched$CHROM)
  outMatrix[, 2] <- as.character(inputStitched$START)
  outMatrix[, 3] <- as.character(inputStitched$STOP)
  outMatrix[, 4] <- as.character(inputStitched$ID)
  
  if (length(score) == nrow(inputStitched)) {
    score <- rank(score, ties.method = "first")
    score <- length(score) - score + 1 
    outMatrix[, 5] <- as.character(score)
  }
  
  trackDescription <- paste(trackDescription, "\nCreated on ", format(Sys.time(), "%b %d %Y"), sep = "")
  trackDescription <- gsub("\n", "\t", trackDescription)
  tName <- gsub(" ", "_", trackName)
  cat('track name="', tName, '" description="', trackDescription, '" itemRGB=On color=', baseColor, "\n", 
      sep = "", file = outputFile)
  write.table(file = outputFile, outMatrix, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  if (splitSuper) {
    cat("\ntrack name=\"Super_", tName, '" description="Super ', trackDescription, '" itemRGB=On color=', superColor, "\n", 
        sep = "", file = outputFile, append = TRUE)
    write.table(file = outputFile, outMatrix[superRows, ], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}
convert_stitched_to_gateway_bed <- function(inputStitched, outputFileRoot, splitSuper = TRUE, 
                                            score = c(), superRows = c()) {
  outMatrix <- matrix(data = "", ncol = 6, nrow = nrow(inputStitched))
  
  outMatrix[, 1] <- as.character(inputStitched$CHROM)
  outMatrix[, 2] <- as.character(inputStitched$START)
  outMatrix[, 3] <- as.character(inputStitched$STOP)
  outMatrix[, 4] <- as.character(inputStitched$ID)
  
  if (length(score) == nrow(inputStitched)) {
    score <- rank(score, ties.method = "first")
    score <- length(score) - score + 1
    outMatrix[, 5] <- as.character(score)
  }
  outMatrix[, 6] <- as.character(rep('.', nrow(outMatrix)))
  outputFile1 <- paste(outputFileRoot, '_Gateway_sequences.bed', sep = '')
  write.table(file = outputFile1, outMatrix, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  if (splitSuper) {
    outputFile2 <- paste(outputFileRoot, '_Gateway_m6Agetsequences.bed', sep = '')
    write.table(file = outputFile2, outMatrix[superRows, ], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}
write6mAget_table <- function(sequence6mAget, description, outputFile, additionalData = NA) {
  
  description <- paste("#", description, "\nCreated on ", format(Sys.time(), "%b %d %Y"), sep = "")
  description <- gsub("\n", "\n#", description)
  cat(description, "\n", file = outputFile)
  
  if (!is.null(additionalData) && is.matrix(additionalData)) {
    if (nrow(additionalData) != nrow(sequence6maget)) {
      stop("The number of rows of the additional data does not match the number of rows of the sequence data. Please check the input!")
    } else {
      sequence6mAget <- cbind(sequence6mAget, additionalData)
      if ("6mARank" %in% colnames(additionalData)) {
        sequence6mAget <- sequence6mAget [order(sequence6mAget $Rank), ]
      }
    }
  } else if (!is.null(additionalData)) {
    warning("If the additional data is not a matrix or is empty, the processing of additional data will be skipped.")
  }
  
  write.table(file = outputFile, sequence6mAget, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
}
outFolder <- "./"
signalsumFile <- ""
sampleName <- ""
stitched_regions <- read.delim(file = signalsumFile, sep = "\t", header = TRUE)
# Column name: ID, CHROM, START, STOP, SIGNAL
signal_vector <- as.numeric(stitched_regions$SIGNAL)
signal_vector[signal_vector < 0] <- 0
cutoff_options <- calculate_cutoff(signal_vector, drawPlot = FALSE)
m6AgetsequenceRows <- which(signal_vector > cutoff_options$absolute)
typical6mA_regions <- setdiff(1:nrow(stitched_regions), m6AgetsequenceRows)
m6Aregions_Description <- paste(
  sampleName, " m6A_get\nCreated from ", signalsumFile, 
  "\nUsing cutoff of ", cutoff_options$absolute, " for m6A_get_sequence", sep = ""
)
plotFileName <- paste(outFolder, sampleName, '_Plot_points.png', sep = '')
png(filename = plotFileName, height = 600, width = 600)
signalOrder <- order(signal_vector, decreasing = TRUE)
plot(length(signal_vector):1, signal_vector[signalOrder], col = 'red', 
     xlab = '6mA-regions', ylab = 'Signal', pch = 19, cex = 2)
abline(h = cutoff_options$absolute, col = 'grey', lty = 2)
abline(v = length(signal_vector) - length(m6AgetsequenceRows), col = 'grey', lty = 2)
lines(length(signal_vector):1, signal_vector[signalOrder], lwd = 4, col = 'red')
text(0, 0.8 * max(signal_vector), paste(
  ' Cutoff used: ', cutoff_options$absolute, '\n',
  '6mA-get-sequence identified: ', length(m6AgetsequenceRows)
), pos = 4)
dev.off()
