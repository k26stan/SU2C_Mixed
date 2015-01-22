## SU2C Data: Logistic Fit, permuting for significance ##
## Mixed Effects Models vs Gene Expression ##
## Kristopher Standish ##
## January 7, 2014 ##

#################################################################
## POTENTIAL DRUGS ##############################################
#################################################################

## Old Version of what to plot
# TO_PLOT <- c(4,11,17,18,19,21,23,33,41,47,49,64,66,69,73,79,96,100,101,112,119,122,125,126)

## Potential Drugs of Interest
DRUGS.pref <- c("Trametinib","MEK162","Palbociclib","MLN0128","GSK2141795")
DRUGS.pref2 <- c("Trametinib","MEK-162","PD325901","Palbociclib","MLN0128","INK_128","OSI-027","everolimus","sirilimus","temsirilimus","GSK2141795")
DRUGS.ccle <- c("Topotecan","Nilotinib","Lapatinib","Irinotecan","Erlotinib","Sorafenib")
DRUGS.prev <- c(4,11,17,18,19,21,23,33,41,47,49,64,66,69,73,79,96,100,101,112,119,122,125,126)
DRUGS.fit <- c(11,18,33,47,69,122,126)
DRUGS.gen <- c("Cabozantinib","Dacomitinib","Etoposide","MLN2480","Palbociclib","Sunitinib","Vorinostat")

#################################################################
## GET ORGANIZED ################################################
#################################################################

## Set Date/Identifier
DATE <- "20150109"

## Load Packages
library(lcmm)
library(lattice)
library(gplots)
library(ggplot2)
library(nlme)
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
## LOAD GENE EXPRESSION DATA ####################################
#################################################################

TOP <- GTOP <- list()

## Load Top Genes and Expression Data

 # Top Gene Names
TOP$dox <- read.table(paste(PathToTop2,"doxorubicin_top50_Genes_allcelllines.txt",sep=""),sep="\t",header=T)
names(TOP$dox) <- c("Probe","Gene","P")
TOP$iri <- read.table(paste(PathToTop3,"Irinotecan_top50_Genes_allcelllines.txt",sep=""),sep="\t",header=T)
names(TOP$iri) <- c("Probe","Gene","P")
TOP$top <- read.table(paste(PathToTop3,"Topotecan_top50_Genes_allcelllines.txt",sep=""),sep="\t",header=T)
names(TOP$top) <- c("Probe","Gene","P")

 # Top Expression Values
GTOP$dox.1 <- read.table(paste(PathToTop2,"doxorubicin_top50_GEX_allcelllines.txt",sep=""),sep="\t",header=T)
names(GTOP$dox.1)[1] <- c("Probe")
GTOP$iri.1 <- read.table(paste(PathToTop3,"Irinotecan_top50_GEX_allcelllines.txt",sep=""),sep="\t",header=T)
names(GTOP$iri.1)[1] <- c("Probe")
GTOP$top.1 <- read.table(paste(PathToTop3,"Topotecan_top50_GEX_allcelllines.txt",sep=""),sep="\t",header=T)
names(GTOP$top.1)[1] <- c("Probe")

#################################################################
## ORGANIZE GENE EXPRESSION DATA ################################
#################################################################

## Specify Drug Names
ALL_DRUGS <- colnames(DAT.1)[7:134]
BTOP <- KEY <- GEX <- KEYS <- list()
 # and Shorthand
DRUG_SHORT_HAND <- c("dox","iri","top")
names(DRUG_SHORT_HAND) <- c( "DoxorubicinHCl","IrinotecanHCl","TopotecanHCl")

## Loop through Drugs of Interest
for ( drug in DRUG_SHORT_HAND ) {
	## Create Gene-Probe Key
	KEY[[drug]] <- TOP[[drug]][which(!duplicated(TOP[[drug]]$Probe)),c("Gene","Probe")]
	KEY[[drug]] <- data.frame(KEY[[drug]],TAG=paste(KEY[[drug]]$Gene,KEY[[drug]]$Probe,sep="_"))
	KEY[[drug]]$TAG <- gsub("---",".",KEY[[drug]]$TAG, fixed=T)
	KEY[[drug]]$TAG <- gsub("-",".",KEY[[drug]]$TAG, fixed=T)
	KEY[[drug]]$TAG <- gsub(" ","",KEY[[drug]]$TAG, fixed=T)
	KEY[[drug]]$TAG <- gsub("///",".",KEY[[drug]]$TAG, fixed=T)	
	
	## Give Expression Values Gene Names
	GTOP[[drug]] <- merge(KEY[[drug]],GTOP[[paste(drug,".1",sep="")]])
	BTOP[[drug]] <- data.frame(GTOP[[drug]][,4:ncol(GTOP[[drug]])])
	for ( which_gene in 1:nrow(BTOP[[drug]]) ) {
		BTOP[[drug]][which_gene,] <- as.numeric( as.numeric(BTOP[[drug]][which_gene,]) > median(as.numeric(BTOP[[drug]][which_gene,])) )
	} # apply(GTOP[,2:ncol(GTOP)],2,median) )
	rownames(BTOP[[drug]]) <- GTOP[[drug]]$TAG

	## Compile GEX data into list
	GEX[[drug]] <- t( BTOP[[drug]] )
	KEYS[[drug]] <- data.frame(GTOP[[drug]][,1:3],PAREN=paste(GTOP[[drug]][,2]," (",GTOP[[drug]][,1],")",sep=""))
}

#################################################################
## USE NLME TO FIT LOGISTIC MODELS ##############################
#################################################################

## Compile a few values
DOSES <- unique(DAT.2$Dose_Value)
CELL_LINES <- as.character(unique(DAT.2$Cell))

## Make Function to Plot this Shiz
GEXRUN <- function(drug_name,start,stop) {
	## Specify Drug Name
	DOI <- drug_name
	DOI.short <- DRUG_SHORT_HAND[DOI]
	DOI2 <- DAT.2[,DOI]
	print(paste("%%%%%%%% Drug =",DOI,"%%%%%%%%"))

	## Set up New Data Frame w/ Only Drug of Interest
	COLS_PULL <- c(names(DAT.2)[1:6],"Dose","iline","Batch",DOI)
	DAT.3a <- DAT.2[,COLS_PULL]
	names(DAT.3a)[ncol(DAT.3a)] <- "Resp"

	## Put together new DAT.3 data.frame
	TEMP_GEX <- GEX[[DOI.short]][match(DAT.3a$Cell,rownames(GEX[[DOI.short]])),]
	DAT.3b <- data.frame(DAT.3a, TEMP_GEX )
	# Getting Pos/Neg controls for Permutation testing
	# for ( col in 11:60 ) { print(paste( "##",col-10,"-", formatC(summary(lm(DAT.3b$Resp ~ DAT.3b[,col]))$coefficients[2,4],format="e",digits=1) )) }
	# (Dox+Ctrl==26; Dox-Ctrl==8)

	## Set Initial Estimates for Nonlinear Fit
	START <- c(Asym=100, xmid=0, scal=-.5)
	start_time <- proc.time()
	if (stop==0) { stop <- ncol(TEMP_GEX) }
	for ( probe_ind in start:stop ) {
		## Specify Probe to Check Out
		PROBE <- colnames(TEMP_GEX)[probe_ind]
		PROBE_PAR <- KEYS[[drug_name]][which(KEYS[[drug_name]]$TAG==PROBE),"PAREN"]
		print(paste("%%%%%%%% Probe",probe_ind,"=",PROBE,"%%%%%%%%"))

		## Make data frame with only 1 Probe
		DAT.3 <- data.frame(DAT.3b[,1:which(colnames(DAT.3b)=="Resp")],DAT.3b[,PROBE])
		names(DAT.3)[ncol(DAT.3)] <- "Probe"
		DAT.3 <- DAT.3[which(!is.na(DAT.3$Probe)),]

		###############################################
		## Random = ~ xmid | iline
		print("%%%%%%%% Random = ~ xmid | iline %%%%%%%%")

		## Run Nonlinear Fit for all cells (using xmid as Random Effect)
		print("Running MOD_ALL")
		MOD_ALL <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | iline, data=DAT.3, start=START, subset=!is.na(Probe) ),T)
		if (is.character(MOD_ALL)==F) {
			EST_ALL <- coef(MOD_ALL)[[1]]
			SUM_ALL <- summary(MOD_ALL)
			FIT_ALL <- SUM_ALL$coefficients
			LL_ALL <- logLik(MOD_ALL)[1]
			PAR_ALL <- attr(logLik(MOD_ALL),"df")
		}else { EST_ALL <- SUM_ALL <- FIT_ALL <- list() ; LL_ALL <- PAR_ALL <- 0  }
		## Run Nonlinear Fit for BRAF Mutant Status
		 # Mutant
		print("Running MOD_Gup")
		MOD_Gup <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | iline, data=DAT.3, start=START, subset=Probe==1 ),T)
		if (is.character(MOD_Gup)==F) {
			EST_Gup <- coef(MOD_Gup)[[1]]
			SUM_Gup <- summary(MOD_Gup)
			FIT_Gup <- SUM_Gup$coefficients
			LL_Gup <- logLik(MOD_Gup)[1]
			PAR_Gup <- attr(logLik(MOD_Gup),"df")
		}else { EST_Gup <- SUM_Gup <- FIT_Gup <- list() ; LL_Gup <- PAR_Gup <- 0  }
		 # Wildtype
		print("Running MOD_Gdn")
		MOD_Gdn <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | iline, data=DAT.3, start=START, subset=Probe==0 ),T)
		if (is.character(MOD_Gdn)==F) {
			EST_Gdn <- coef(MOD_Gdn)[[1]]
			SUM_Gdn <- summary(MOD_Gdn)
			FIT_Gdn <- SUM_Gdn$coefficients
			LL_Gdn <- logLik(MOD_Gdn)[1]
			PAR_Gdn <- attr(logLik(MOD_Gdn),"df")
		}else { EST_Gdn <- SUM_Gdn <- FIT_Gdn <- list() ; LL_Gdn <- PAR_Gdn <- 0  }

		## Calculate Likelihood Ratio for 1 vs 2
		P_VAL_G <- 1
		if ( LL_ALL!=0 & LL_Gup!=0 & LL_Gdn!=0 ) {
			LR_ALL_G <- -2*( LL_ALL - (LL_Gup + LL_Gdn) )
			DF_ALL_G <- PAR_Gup + PAR_Gdn - PAR_ALL # [(LLp+ + LLp-) - LLpt] 
			P_VAL_G <- dchisq(LR_ALL_G, DF_ALL_G)
		}

		###############################################
		## PERMUTE ## Random = ~ xmid | iline
		n_perm <- 100
		LL.perm <- array( , c(n_perm,2) ) ; colnames(LL.perm) <- c("Gup","Gdn")
		PAR.perm <- array( , c(n_perm,2) ) ; colnames(PAR.perm) <- c("Gup","Gdn")
		LR.perm <- P.perm <- DF.perm <- numeric( n_perm )
		print("%%%%%%%% Starting Permutations %%%%%%%%")
		# for ( p in 1:n_perm ) {
		keep_permuting <- 1
		p <- 0
		start_perm <- proc.time()
		while ( keep_permuting==1 ) {
			p <- p + 1
			print(paste( "Perm:",p,"-",round(proc.time()-start_perm,3)[3] ))
			## Create Permuted Response Table (retaining Dosing)
			DAT.3.perm <- DAT.3
			DAT.3.perm[,"Probe"] <- 0
			WHICH_CELLS <- sample( CELL_LINES, length(CELL_LINES)/2, replace=F )
			WHICH_ROWS <- which( DAT.3.perm[,"Cell"] %in% WHICH_CELLS )
			DAT.3.perm[WHICH_ROWS,"Probe"] <- 1
			
			## Re-Run Models ##
			 # GEX Up
			print("Up")
			MOD_Gup.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | iline, data=DAT.3.perm, start=START, subset=Probe==1 ),T)
			if (is.character(MOD_Gup.perm)==F) {
				EST_Gup.perm <- coef(MOD_Gup.perm)[[1]]
				SUM_Gup.perm <- summary(MOD_Gup.perm)
				FIT_Gup.perm <- SUM_Gup.perm$coefficients
				LL_Gup.perm <- logLik(MOD_Gup.perm)[1]
				PAR_Gup.perm <- attr(logLik(MOD_Gup.perm),"df")
			}else { EST_Gup.perm <- SUM_Gup.perm <- FIT_Gup.perm <- list() ; LL_Gup.perm <- PAR_Gup.perm <- 0  }
			 # GEX Down
			print("Down")
			MOD_Gdn.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | iline, data=DAT.3.perm, start=START, subset=Probe==0 ),T)
			if (is.character(MOD_Gdn.perm)==F) {
				EST_Gdn.perm <- coef(MOD_Gdn.perm)[[1]]
				SUM_Gdn.perm <- summary(MOD_Gdn.perm)
				FIT_Gdn.perm <- SUM_Gdn.perm$coefficients
				LL_Gdn.perm <- logLik(MOD_Gdn.perm)[1]
				PAR_Gdn.perm <- attr(logLik(MOD_Gdn.perm),"df")
			}else { EST_Gdn.perm <- SUM_Gdn.perm <- FIT_Gdn.perm <- list() ; LL_Gdn.perm <- PAR_Gdn.perm <- 0  }

			## Compile Permutation Results
			LL.perm[p,] <- c( LL_Gup.perm, LL_Gdn.perm)
			PAR.perm[p,] <- c( PAR_Gup.perm, PAR_Gdn.perm)
			LR.perm[p] <- -2*( LL_ALL - (LL_Gup.perm + LL_Gdn.perm) )
			 # In the LRT, degrees of freedom is equal to the number of additional parameters in the more complex model.
			DF.perm[p] <- sum( PAR.perm[p,] )
			P.perm[p] <- dchisq(LR.perm[p], DF.perm[p])
			## If it's clearly not going to be significant, stop permuting
			h <- length( which(LR.perm > LR_ALL_G ))
			if ( h >= 1 | p == n_perm ) { keep_permuting <- 0 }
		}
		P_PERM <- ( length(which( LR.perm[1:p] > LR_ALL_G )) + 1 ) / (p+1)

		###########################################
		## Plot Fits based on Model Coefficients ##
		PCH_COLS <- c("mediumpurple3","sienna3")
		FIT_COLS <- c("mediumpurple1","sienna1","mediumpurple4","sienna3","black")
		X_VALS <- seq(-2,2,.01)
		png(paste(PathToSave,"PL_1-NLME_iline_Gene_",DOI,"-",PROBE,".png",sep=""), width=2000,height=1000,pointsize=24)
		par(mfrow=c(1,2))
		## Plot Real Data
		print("### Plotting True Data ###")
		print("Plotting Empty Plot and MOD_ALL")
		plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("GEx & Dose-Resp - ",DOI,": ",PROBE_PAR,sep=""))
		abline( h=seq(0,200,20), col="grey50", lty=c(1,2,2,2,2,1,2,2,2,2,1) )
		abline( v=seq(-3,3,1), col="grey50", lty=2 )
		points(DAT.3[,"Dose"],DAT.3[,"Resp"], pch="+", col=PCH_COLS[factor(DAT.3$Probe)] )
		GEX_COLS <- FIT_COLS[factor(DAT.3$Probe[seq(1,540,9)])]
		if (is.character(MOD_ALL)==F) {
			for (i in 1:nrow(EST_ALL) ) {
				Y_VALS <- EST_ALL$Asym[i] / ( 1 + exp( -(X_VALS - EST_ALL$xmid[i]) / EST_ALL$scal[i] ) )	
				points(X_VALS,Y_VALS, type="l", col=GEX_COLS[i] )
			}
			Y_VALS <- FIT_ALL["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_ALL["xmid","Estimate"]) / FIT_ALL["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[5], lwd=6, lty=2)
		}
		print("Plotting Gup")
		if (is.character(MOD_Gup)==F) {
			Y_VALS <- FIT_Gup["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gup["xmid","Estimate"]) / FIT_Gup["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[4], lwd=6, lty=2)
		}
		print("Plotting Gdn")
		if (is.character(MOD_Gdn)==F) {
			Y_VALS <- FIT_Gdn["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gdn["xmid","Estimate"]) / FIT_Gdn["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[3], lwd=6, lty=2)
		}
		if ( LL_ALL!=0 & LL_Gup!=0 & LL_Gdn!=0 ) {
			text(-2,10, pos=4, labels=paste("Permutation-Based P-Value (",p,")",sep=""), col="black")
			text(-2,5, pos=4, labels=round(P_PERM,4) ) # paste("p =",formatC(P_PERM,2,format="e")) )
		}
		text(1.2,115, pos=4, labels="(True Data)", col="black")
		## Plot Randomly Permuted Data
		print("### Plotting Permuted Version of Data ###")
		print("Plotting Empty Plot and MOD_ALL")
		plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("GEx & Dose-Resp - ",DOI,": ",PROBE_PAR,sep=""))
		abline( h=seq(0,100,20), col="grey50", lty=c(1,2,2,2,2,1) )
		abline( v=seq(-3,3,1), col="grey50", lty=2 )
		points(DAT.3.perm[,"Dose"],DAT.3.perm[,"Resp"], pch="+", col=PCH_COLS[factor(DAT.3.perm$Probe)] )
		GEX_COLS <- FIT_COLS[factor(DAT.3.perm$Probe[seq(1,540,9)])]
		if (is.character(MOD_ALL)==F) {
			for (i in 1:nrow(EST_ALL) ) {
				Y_VALS <- EST_ALL$Asym[i] / ( 1 + exp( -(X_VALS - EST_ALL$xmid[i]) / EST_ALL$scal[i] ) )	
				points(X_VALS,Y_VALS, type="l", col=GEX_COLS[i] )
			}
			Y_VALS <- FIT_ALL["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_ALL["xmid","Estimate"]) / FIT_ALL["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[5], lwd=6, lty=2)
		}
		print("Plotting Gup")
		if (is.character(MOD_Gup.perm)==F) {
			Y_VALS <- FIT_Gup.perm["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gup.perm["xmid","Estimate"]) / FIT_Gup.perm["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[4], lwd=6, lty=2)
		}
		print("Plotting Gdn")
		if (is.character(MOD_Gdn.perm)==F) {
			Y_VALS <- FIT_Gdn.perm["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gdn.perm["xmid","Estimate"]) / FIT_Gdn.perm["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[3], lwd=6, lty=2)
		}
		legend(.5,110, legend=c("Ind - GEx Up","Ind - GEx Down","Fit - ALL","Fit - GEx Up","Fit - GEx Down"),lty=c(1,1,2,2,2),col=FIT_COLS[c(2,1,5,4,3)], lwd=c(1,1,6,6,6) )		# if ( LL_ALL!=0 & LL.perm[p,2]!=0 & LL.perm[p,3]!=0 ) { text(-2,20, pos=4, labels=paste("p =",formatC(P_VAL_G,2,format="e")), col="black") }
		text(.4,115, pos=4, labels="(Randomly Permuted Data)", col="black")
		dev.off()

	} # Close for loop of Probes
} # Close GEXRUN function

#################################################################
## RUN IT #######################################################
############# GEXRUN(drug_name,start,stop) ######################

GEXRUN("DoxorubicinHCl",1,0)
# GEXRUN("Dacarbazine",1,0)
# GEXRUN("Sorafenib",1,0)
# GEXRUN("Etoposide",1,0)
# GEXRUN("Vorinostat",1,0)
# GEXRUN("CabozantinibXL184",1,0)
# GEXRUN("DacomitinibPF299804",1,0)
# GEXRUN("MLN2480",1,0)
# GEXRUN("Sunitinib",1,0)
# GEXRUN("PalbociclibPD0332991Isethionate",1,0)

## Look into:
#1) Trametinib - 103
#1b) MEK162 (????) - 108
#2) Palbociclib - 106
#3) MLN0128 - 105
#4) GSK2141795 - ??????



