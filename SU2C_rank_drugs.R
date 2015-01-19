## SU2C Data: Ranks Drugs based on Response Data ##
## Kristopher Standish ##
## January 7, 2014 ##

## CCLE logistic fits (Kuan)
## 4 parameter fits
## Fit vs GEX
## Rank Drugs by quality fit
 # IQR
 # LogLik
## 

#################################################################
## GET ORGANIZED ################################################
#################################################################

## Set Date/Identifier
DATE <- "20150116"

## Load Packages
library(lcmm)
library(lattice)
library(gplots)
library(ggplot2)
# library(nlme)
library(lme4)
library(xlsx)

## Load Time-Series Data
# Mac Paths
PathToData <- "/Users/kstandis/Data/SU2C/"
PathToTop <- "/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/From_Kuan/"
PathToTop2 <- "/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/From_Kuan_2/"
PathToTop3 <- "/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/From_Kuan_3/"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/",DATE,"_Plots/",sep="")

#################################################################
## LOAD DOSE-RESPONSE DATA ######################################
#################################################################

## Load Previously Compiled Dose-Response Data, missing several cell lines
DAT.1 <- read.table(paste(PathToData,"alldrugs_subset_probeset.txt",sep=""),sep="\t",header=T)
colnames(DAT.1)[c(1,6)] <- c("Cell","Dose_Value")
DAT.1 <- data.frame(DAT.1,Dose=log10(DAT.1$Dose_Value))
colnames(DAT.1) <- gsub(".","",colnames(DAT.1),fixed=T)

## Load Updated Compiled Dose-Response Data, with all cell lines
DAT.1b <- read.table(paste(PathToData,"SU2C_alldrugs_allcellines_9pt.txt",sep=""),sep="\t",header=T)
colnames(DAT.1b)[c(1,6)] <- c("Cell","Dose_Value")
DAT.1b <- data.frame(DAT.1b,Dose=log10(DAT.1b$Dose_Value))
colnames(DAT.1b) <- gsub(".","",colnames(DAT.1b),fixed=T)

## Load New Tables for Cell Lines
FILE_LIST <- as.character( read.table(paste(PathToData,"20141201_Filenames.txt",sep=""), header=F, sep="\t" )[,1] )
LNS.dat <- LNS.dos <- list()
start <- proc.time()
for ( i in 1:length(FILE_LIST) ) {
	file_name <- FILE_LIST[i]
	short_name <- gsub("-","", strsplit( file_name, " " )[[1]][1])
	LNS.dat[[short_name]] <- read.xlsx( paste(PathToData,"20141201_Cell_Lines/",file_name,sep=""), sheetIndex=1, colIndex=13:40, rowIndex=11:138, header=F, keepFormulas=F)
	LNS.dos[[short_name]] <- read.xlsx( paste(PathToData,"20141201_Cell_Lines/",file_name,sep=""), sheetIndex=1, colIndex=14:40, rowIndex=9, header=F, keepFormulas=F)
	print(paste( "Done w/ ",i,"of",length(FILE_LIST),"-",round((proc.time()-start)[3],1) ))
}

## Get list of Old/New Cell Lines
OLD_LINES <- as.character( unique( DAT.1$Cell ) )
NEW_LINES <- setdiff( as.character(unique(DAT.1b$Cell)), OLD_LINES )

#################################################################
## FORMAT NEW DOSE-RESPONSE TABLES ##############################
#################################################################

## Specify Dose, Rep, etc...
# DOSE <- unique( DAT.1$Dose_Value )[c( rep(c(1,2,3),3), rep(c(4,5,6),3), rep(c(7,8,9),3) )]
REP <- rep( rep( c(1,2,3), rep(3,3) ), 3)

## Compile/Reformat New Data
start <- proc.time()
for ( i in 1:length(LNS.dat) ) {
	CELL <- names(LNS.dat)[i]
	LEN <- length( LNS.dos[[i]] )
	MTA.temp <- data.frame( Cell=rep(CELL,LEN), BRAF=rep(0,LEN), NRAS=rep(0,LEN), KRAS=rep(0,LEN), rep=REP, Dose_Value=c(LNS.dos[[i]],recursive=T) )
	DAT.temp <- t( LNS.dat[[i]][,2:ncol(LNS.dat[[i]])] )
	colnames(DAT.temp) <- LNS.dat[[i]][,1]
	if ( i == 1 ) {
		LNS.1 <- data.frame( MTA.temp, DAT.temp )
	}else{
		COMP.temp <- data.frame( MTA.temp, DAT.temp )
		LNS.1 <- rbind( LNS.1, COMP.temp )
	}
	print(paste( "Done w/",i,"of",length(LNS.dat),"-",round((proc.time()-start)[3],1) ))
}
colnames(LNS.1) <- gsub(".","",colnames(LNS.1), fixed=T)
LNS.2 <- data.frame( LNS.1, random_probeset=rep(0,nrow(LNS.1)), Dose=log10(LNS.1$Dose_Value) )
identical( colnames(LNS.2), colnames(DAT.1) )

## Make DAT.comp from Original Data + Manually Compiled Data
DAT.comp.1 <- rbind( DAT.1, LNS.2 )
 # Filter out superfluous columns
DAT.comp <- DAT.1b[ , which( colnames(DAT.1b) %in% colnames(DAT.1) ) ]
 # Make Indicator for different cell lines
iline1 <- paste(DAT.comp$Cell,DAT.comp$rep,sep="_")
iline <- as.numeric(factor(iline1))
 # Include "iline" and "Batch" in Data Frame
DAT.2 <- data.frame( DAT.comp, iline, Batch=rep(2,nrow(DAT.comp)) )
BATCH_1 <- which( DAT.comp$Cell %in% OLD_LINES )
DAT.2$Batch[BATCH_1] <- 1

#################################################################
## FCT: RANK DRUGS FOR QUALITY OF DATA ##########################
#################################################################

## Compile a few values
DOSES <- unique(DAT.2$Dose_Value)
CELL_LINES <- as.character(unique(DAT.2$Cell))
DRUG_LIST <- names(DAT.2)[7:134]
drug_list <- DRUG_LIST[ c(11,18,33,47,69,122,126) ]

## Make Function to Rank Drugs Quickly
RANK.fast <- function(drug_list) {
	COMP.iqr <- array( , c(length(drug_list),8) )
	colnames(COMP.iqr) <- c( names(summary(DAT.2[,7])), "IQR" )
	COMP.lin <- array( , c(length(drug_list),3) )
	colnames(COMP.lin) <- c( "COEF","P","LL" )
	COMP.log <- array( , c(length(drug_list),6) )
	colnames(COMP.log) <- c( "Asym_1","Asym_2","Diff","LL","Est_xmid","Est_scal" )
	rownames(COMP.iqr) <- rownames(COMP.lin) <- rownames(COMP.log) <- drug_list
	d <- 0
	for ( DOI in drug_list ) {
		d <- d + 1
		# DOI <- "DoxorubicinHCl"
		# DOI <- "Cladribine"
		## Specify Drug Name
		DOI2 <- DAT.2[,DOI]
		print(paste("%%%%%%%% Drug =",DOI,"%%%%%%%%"))

		## Set up New Data Frame w/ Only Drug of Interest
		COLS_PULL <- c(names(DAT.2)[1:6],"Dose","Cell","Batch",DOI)
		DAT.3a <- DAT.2[,COLS_PULL]
		names(DAT.3a)[ncol(DAT.3a)] <- "Resp"
		DAT.3 <- DAT.3a

		###############################################
		## Rank based on IQR
		SUMMARY <- summary( DAT.3[,"Resp"] )
		COMP.iqr[d,names(SUMMARY)] <- SUMMARY
		COMP.iqr[d,"IQR"] <- IQR(DAT.3[,"Resp"], na.rm=T)

		###############################################
		## Rank based on Linear Fit
		MOD <- lm( Resp ~ Dose, data=DAT.3 )
		COEF <- coef(MOD)["Dose"]
		P <- summary(MOD)$coefficients["Dose","Pr(>|t|)"]
		LL <- logLik(MOD)
		COMP.lin[d,] <- c( COEF, P, LL )

		# ###############################################
		# ## Rank based on Logistic Fit
		START <- c(Asym_1=0, Asym_2=100, xmid=0, scal=.5)
		MOD <- try( nls( Resp ~ Asym_2 + (Asym_1-Asym_2) / (1+exp((xmid-Dose)/scal)), data=DAT.3, start=START ) ,T)
		if (is.character(MOD)==F) {
			ASYM_1 <- coef(MOD)["Asym_1"]
			ASYM_2 <- coef(MOD)["Asym_2"]
			XMID <- coef(MOD)["xmid"]
			SCAL <- coef(MOD)["scal"]
			# SCAL <- coef(MOD)["Asym"]
			# P <- summary(MOD)$coefficients["Asym","Pr(>|t|)"]
			LL <- logLik(MOD)
			COMP.log[d,] <- c( ASYM_1, ASYM_2, ASYM_2-ASYM_1, LL, XMID, SCAL )
		}

	} # Close Drug Loop

	## Compile Outputs
	COMPILE <- list( COMP.iqr, COMP.lin, COMP.log )
	names(COMPILE) <- c("IQR","LIN","LOG")
	return(COMPILE)
} # Close "RANK.fast" function
# RANKOUT.1 <- RANK.fast( DRUG_LIST )

#################################################################
## FCT: RANK DRUGS FOR QUALITY OF MIXED MODEL FIT ###############
#################################################################

## Make Function to Rank Drugs Slowly (based on mixed models)
RANK.slow <- function(drug_list) {
	COMP.mem <- array( , c(length(drug_list),5) )
	colnames(COMP.mem) <- c( "EST_Asym","EST_xmid","EST_scal","LL","Conv" )
	rownames(COMP.mem) <- drug_list
	d <- 0
	for ( DOI in drug_list ) {
		d <- d + 1
		# DOI <- "DoxorubicinHCl"
		# DOI <- "Cladribine"
		## Specify Drug Name
		DOI2 <- DAT.2[,DOI]
		print(paste("%%%%%%%% Drug =",DOI,"%%%%%%%%"))

		## Set up New Data Frame w/ Only Drug of Interest
		COLS_PULL <- c(names(DAT.2)[1:6],"Dose","Cell","Batch",DOI)
		DAT.3a <- DAT.2[,COLS_PULL]
		names(DAT.3a)[ncol(DAT.3a)] <- "Resp"
		DAT.3 <- DAT.3a

		## Set Initial Estimates for Nonlinear Fit
		START <- c(Asym=100, xmid=0, scal=-.5)
		start_time <- proc.time()

		###############################################
		## Random = ~ xmid | Cell
		print("%%%%%%%% Random = ~ xmid | Cell %%%%%%%%")
		 # Run Nonlinear Fit for all cells (using xmid as Random Effect)
		MOD_ALL <- WARN <- ""
		MOD_ALL <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START ),T)
		if (is.character(MOD_ALL)==F) {
			EST_ALL <- coef(MOD_ALL)[[1]]
			SUM_ALL <- summary(MOD_ALL)
			FIT_ALL <- SUM_ALL$coefficients
			LL_ALL <- logLik(MOD_ALL)[1]
			PAR_ALL <- attr(logLik(MOD_ALL),"df")
			# Compile Estimates
			COMP.mem[d,1:4] <- c( FIT_ALL[,1], LL_ALL )
		}else {
			EST_ALL <- SUM_ALL <- FIT_ALL <- list()
			LL_ALL <- PAR_ALL <- 0
		}
		 # Check for convergence
		WARN <- warnings()
		if ( any( grepl("failure to converge",WARN)) ) {
			COMP.mem[,"Conv"] <- 0
		}else{ COMP.mem[,"Conv"] <- 1 }

		## DONE ##
	} # Close Drug Loop
} # Close "RANK.slow" Function
# RANKOUT.2 <- RANK.slow( DRUG_LIST )

#################################################################
## RANK DRUGS BASED ON FUNCTIONS ################################
#################################################################

## Do Quick Ranking for ALL Drugs
WHICH_LIST <- DRUG_LIST
RANKOUT.1 <- RANK.fast( WHICH_LIST )
n.DRUGS <- length( WHICH_LIST )

## Pull out data.frames
 # IQR
IQR.1 <- RANKOUT.1$IQR
 # Linear Model
LIN.1 <- RANKOUT.1$LIN
LIN.1[,2] <- -log10(LIN.1[,2])
 # Logistic Model
LOG.1 <- RANKOUT.1$LOG
LOG.1[which(is.na(LOG.1[,"Diff"])),"Diff"] <- 0

## Compile Relevant Data
COMP.rel <- data.frame(IQR.1[,c("1st Qu.","IQR")],LIN.1[,c("COEF","P")],ASYM=LOG.1[,"Diff"])
 # Plot relevant values against one another
pairs( COMP.rel, col="tomato2", pch="+" )


## Plot ranks
COLS <- c("firebrick2","chocolate2","gold2","chartreuse1","deepskyblue1","mediumpurple3")
XLIM <- c(0,n.DRUGS)

## Log Asym
YLIM <- c(0,max(LOG.1[,1]))
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM )
points( 1:n.DRUGS, sort(LOG.1[,1],decreasing=T) )

## Log Asym
YLIM <- c(0,max(LOG.1[,1]))
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM )
points( 1:n.DRUGS, sort(LOG.1[,1],decreasing=T) )

for ( i in 1:ncol(LOG.1) ) {
	SCALED <- LOG.1[,i] / max(LOG.1[,i],na.rm=T)
	points( 1:n.DRUGS, sort(SCALED,decreasing=T) )
}
	








#################################################################
## END OF DOC ###################################################
#################################################################
