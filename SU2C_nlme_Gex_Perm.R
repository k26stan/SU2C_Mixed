## SU2C Data: Logistic Fit, permuting for significance ##
## Mixed Effects Models vs Gene Expression ##
## 3 parameter fit
## Kristopher Standish ##
## January 7, 2014 ##

#################################################################
## GET ORGANIZED ################################################
#################################################################

## Set Date/Identifier
DATE <- "20150121"

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
PathToGEX <- "/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/From_Kuan_3/su2c_expression_all_log.txt"
PathToGenes <- "/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/From_Kuan_Drug_Genes/"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/Cancer/SU2C/Mixed_Effects/",DATE,"_Plots/",sep="")

#################################################################
## LOAD DRUG/GENE DATA ##########################################
#################################################################

## Load List of Drug Names
DRUG_FILES <- as.character( read.table( paste(PathToGenes,"DRUG_NAMES.txt",sep=""), header=F )[,1] )

## Load List of Genes/Probe IDs
PROBE_KEY <- read.table( paste(PathToGenes,"SU2C_Target_GenesProbesetIDs.txt",sep=""), header=T, colClass=rep("character",2), sep="\t" )
colnames(PROBE_KEY) <- c("Probe","Gene")

## Load files w/ Relevant Genes for Drugs
GENE_LIST <- list()
for ( d in 1:length(DRUG_FILES) ) {
	drug <- DRUG_FILES[d]
	GENE_LIST[[drug]] <- read.table( paste(PathToGenes,drug,".txt",sep=""), header=T, sep="\t" )
}

## Load GEX table
GEX_TAB <- read.table( PathToGEX, sep="\t",header=T )

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
## ORGANIZE GENE/DRUG/GEX DATA ##################################
#################################################################

## Pull out Probes associated w/ each Gene
PROBE_LIST <- list()
for ( d in 1:length(DRUG_FILES) ) {
	drug <- DRUG_FILES[d]
	PROBE_LIST[[drug]] <- array( , c(0,2) )
	colnames(PROBE_LIST[[drug]]) <- c("Probe","Gene")
	for ( g in 1:nrow(GENE_LIST[[drug]]) ) {
		gene <- as.character( GENE_LIST[[drug]][g,1] )
		gene <- gsub( " ","", gene )
		PROBE_LIST[[drug]] <- rbind( PROBE_LIST[[drug]], PROBE_KEY[grep(gene,PROBE_KEY[,2]),] )
	}
}

## Pull out Expression values for each relevant Probe
GEX_LIST <- list()
for ( d in 1:length(DRUG_FILES) ) {
	drug <- DRUG_FILES[d]
	GEX_LIST[[drug]] <- merge( PROBE_LIST[[drug]], GEX, by.x="Probe",by.y="Probesets" )
}

#################################################################
## ORGANIZE GENE EXPRESSION DATA ################################
#################################################################

## Specify Drug Names
ALL_DRUGS <- colnames(DAT.1)[7:134]
DRUG_KEY <- array( ,c(length(ALL_DRUGS),2) )
colnames(DRUG_KEY) <- c("DR","GEX")
DRUG_KEY[,"DR"] <- ALL_DRUGS
for ( d in 1:length(DRUG_FILES) ) {
	drug <- DRUG_FILES[d]
	DRUG_KEY[ grep(drug,DRUG_KEY[,"DR"]), "GEX"] <- drug
}
DRUG_KEY[ grep("Bosut",DRUG_KEY[,"DR"]), "GEX" ] <- "Bosulif"
DRUG_KEY[ grep("Caboz",DRUG_KEY[,"DR"]), "GEX" ] <- "Bosulif"
DRUG_KEY[ grep("Bosut",DRUG_KEY[,"DR"]), "GEX" ] <- "Cometriq"
DRUG_KEY[ grep("Axit",DRUG_KEY[,"DR"]), "GEX" ] <- "Inlyta"
DRUG_KEY[ grep("Dacomit",DRUG_KEY[,"DR"]), "GEX" ] <- "PF00299804"
DRUG_KEY[ grep("Crizo",DRUG_KEY[,"DR"]), "GEX" ] <- "Xalkori"
DRUG_KEY[ grep("Pacli",DRUG_KEY[,"DR"]), "GEX" ] <- "paclitaxel"
DRUG_KEY[ grep("Sunit",DRUG_KEY[,"DR"]), "GEX" ] <- "sutent"
 # Get rid of Drugs w/ no relevant Genes data
DRUG_KEY.2 <- DRUG_KEY[ which(!is.na(DRUG_KEY[,"GEX"])), ]

#################################################################
## USE NLME TO FIT LOGISTIC MODELS ##############################
#################################################################

## Compile a few values
DOSES <- unique(DAT.2$Dose_Value)
CELL_LINES <- as.character(unique(DAT.2$Cell))

## Make Function to Plot this Shiz
GEXRUN <- function(drug_name,which_probes) {
	# DOI <- "DoxorubicinHCl"
	# DOI <- "Cladribine"
	## Specify Drug Name
	DOI <- drug_name
	DOI.gex <- DRUG_KEY[ which(DRUG_KEY[,"DR"]==DOI), "GEX" ]
	print(paste("%%%%%%%% Drug =",DOI,"%%%%%%%%"))

	## Set up New Data Frame w/ Only Drug of Interest
	COLS_PULL <- c(names(DAT.2)[1:6],"Dose","iline","Batch",DOI)
	DAT.3a <- DAT.2[,COLS_PULL]
	names(DAT.3a)[ncol(DAT.3a)] <- "Resp"

	## Include all GEX data in data frame
	TEMP_GEX <- t( GEX_LIST[[DOI.gex]][,3:ncol(GEX_LIST[[DOI.gex]])] )
	TEMP_KEY <- GEX_LIST[[DOI.gex]][,1:2]
	colnames(TEMP_GEX) <- TEMP_KEY[,"Probe"]
	DAT.3b <- merge( DAT.3a, TEMP_GEX, by.x="Cell",by.y="row.names" )

	## Set Initial Estimates for Nonlinear Fit
	COMPILE.LL <- COMPILE.LR <- COMPILE.LL.perm <- COMPILE.LR.perm <- list()
	START <- c(Asym=100, xmid=0, scal=-.5)
	start_time <- proc.time()
	for ( probe_ind in which_probes ) {
		## Specify Probe to Check Out
		PROBE <- TEMP_KEY[probe_ind,"Probe"]
		PROBE.gene <- TEMP_KEY[probe_ind,"Gene"]
		PROBE.gene <- gsub(" ","",PROBE.gene)
		PROBE.gene <- gsub("///","_",PROBE.gene)

		PROBE.print <- paste( PROBE, " (",PROBE.gene,")",sep="" )

		print(paste("%%%%%%%% Probe",probe_ind,"=",PROBE,"%%%%%%%%"))

		## Make data frame with only 1 Probe
		DAT.3 <- data.frame(DAT.3b[,1:which(colnames(DAT.3b)=="Resp")],DAT.3b[,PROBE])
		names(DAT.3)[ncol(DAT.3)] <- "Probe"
		DAT.3 <- DAT.3[which(!is.na(DAT.3$Probe)),]
		 # Stratify Probe based on Median
		DAT.3[,"Probe"] <- c(0,1)[ factor(DAT.3[,"Probe"] >= median(DAT.3[,"Probe"]) ) ]

		###############################################
		## Random = ~ xmid | Cell
		print("%%%%%%%% Random = ~ xmid | Cell %%%%%%%%")

		## Run Nonlinear Fit for all cells (using xmid as Random Effect)
		print("Running MOD_ALL")
		MOD_ALL <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=!is.na(Probe) ),T)
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
		MOD_Gup <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=Probe==1 ),T)
		if (is.character(MOD_Gup)==F) {
			EST_Gup <- coef(MOD_Gup)[[1]]
			SUM_Gup <- summary(MOD_Gup)
			FIT_Gup <- SUM_Gup$coefficients
			LL_Gup <- logLik(MOD_Gup)[1]
			PAR_Gup <- attr(logLik(MOD_Gup),"df")
		}else { EST_Gup <- SUM_Gup <- FIT_Gup <- list() ; LL_Gup <- PAR_Gup <- 0  }
		 # Wildtype
		print("Running MOD_Gdn")
		MOD_Gdn <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=Probe==0 ),T)
		if (is.character(MOD_Gdn)==F) {
			EST_Gdn <- coef(MOD_Gdn)[[1]]
			SUM_Gdn <- summary(MOD_Gdn)
			FIT_Gdn <- SUM_Gdn$coefficients
			LL_Gdn <- logLik(MOD_Gdn)[1]
			PAR_Gdn <- attr(logLik(MOD_Gdn),"df")
		}else { EST_Gdn <- SUM_Gdn <- FIT_Gdn <- list() ; LL_Gdn <- PAR_Gdn <- 0  }

		## Calculate Likelihood Ratio for 1 vs 2
		if ( LL_ALL!=0 & LL_Gup!=0 & LL_Gdn!=0 ) {
			LR_ALL_G <- -2*( LL_ALL - (LL_Gup + LL_Gdn) )
		}

		###############################################
		## PERMUTE ## Random = ~ xmid | Cell
		n_perm <- 100
		h.lim <- 2
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
			MOD_Gup.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3.perm, start=START, subset=Probe==1 ),T)
			if (is.character(MOD_Gup.perm)==F) {
				EST_Gup.perm <- coef(MOD_Gup.perm)[[1]]
				SUM_Gup.perm <- summary(MOD_Gup.perm)
				FIT_Gup.perm <- SUM_Gup.perm$coefficients
				LL_Gup.perm <- logLik(MOD_Gup.perm)[1]
				PAR_Gup.perm <- attr(logLik(MOD_Gup.perm),"df")
			}else { EST_Gup.perm <- SUM_Gup.perm <- FIT_Gup.perm <- list() ; LL_Gup.perm <- PAR_Gup.perm <- 0  }
			 # GEX Down
			print("Down")
			MOD_Gdn.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3.perm, start=START, subset=Probe==0 ),T)
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
			if ( h >= h.lim | p == n_perm ) { keep_permuting <- 0 }
		}
		P_PERM <- ( length(which( LR.perm[1:p] > LR_ALL_G )) + 1 ) / (p+1)

		###########################################
		## Plot Fits based on Model Coefficients ##
		# PCH_COLS <- c("mediumpurple3","sienna3")
		PCH_COLS <- c("slateblue1","springgreen1")
		# FIT_COLS <- c("mediumpurple1","sienna1","mediumpurple4","sienna3","black")
		FIT_COLS <- c("slateblue1","springgreen1","slateblue3","springgreen3","black")
		X_VALS <- seq(-2,2,.01)
		png(paste(PathToSave,"PL_1-NLME_Cell_GEX_",DOI,"-",PROBE.gene,"_",PROBE,".png",sep=""), width=1250,height=1050,pointsize=24)
		# par(mfrow=c(1,2))
		## Plot Real Data
		print("### Plotting True Data ###")
		print("Plotting Empty Plot and MOD_ALL")
		plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("GEx & Dose-Resp - ",DOI,": ",PROBE.print,sep=""))
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
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[5], lwd=8, lty=2)
		}
		print("Plotting Gup")
		if (is.character(MOD_Gup)==F) {
			Y_VALS <- FIT_Gup["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gup["xmid","Estimate"]) / FIT_Gup["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[4], lwd=8, lty=2)
		}
		print("Plotting Gdn")
		if (is.character(MOD_Gdn)==F) {
			Y_VALS <- FIT_Gdn["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gdn["xmid","Estimate"]) / FIT_Gdn["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS[3], lwd=8, lty=2)
		}
		if ( LL_ALL!=0 & LL_Gup!=0 & LL_Gdn!=0 ) {
			text(-2,10, pos=4, labels=paste("Permutation-Based P-Value (",p,")",sep=""), col="black")
			text(-2,5, pos=4, labels=round(P_PERM,4) ) # paste("p =",formatC(P_PERM,2,format="e")) )
		}
		## Legend
		if ( length(which(DAT.3[,"Resp"]>60)) / nrow(DAT.3) > .5 ) {
			LEG.coords <- c(-2,50)
		}else{ LEG.coords <- c(.75,118) }
		legend( LEG.coords[1],LEG.coords[2], legend=c("Ind - GEx Up","Ind - GEx Down","Fit - ALL","Fit - GEx Up","Fit - GEx Down"),lty=c(1,1,2,2,2),col=FIT_COLS[c(2,1,5,4,3)], lwd=c(1,1,6,6,6) )
		# text(1.2,115, pos=4, labels="(True Data)", col="black")
		dev.off()
		# ## Plot Randomly Permuted Data
		# print("### Plotting Permuted Version of Data ###")
		# print("Plotting Empty Plot and MOD_ALL")
		# plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("GEx & Dose-Resp - ",DOI,": ",PROBE_PAR,sep=""))
		# abline( h=seq(0,100,20), col="grey50", lty=c(1,2,2,2,2,1) )
		# abline( v=seq(-3,3,1), col="grey50", lty=2 )
		# points(DAT.3.perm[,"Dose"],DAT.3.perm[,"Resp"], pch="+", col=PCH_COLS[factor(DAT.3.perm$Probe)] )
		# GEX_COLS <- FIT_COLS[factor(DAT.3.perm$Probe[seq(1,540,9)])]
		# if (is.character(MOD_ALL)==F) {
		# 	for (i in 1:nrow(EST_ALL) ) {
		# 		Y_VALS <- EST_ALL$Asym[i] / ( 1 + exp( -(X_VALS - EST_ALL$xmid[i]) / EST_ALL$scal[i] ) )	
		# 		points(X_VALS,Y_VALS, type="l", col=GEX_COLS[i] )
		# 	}
		# 	Y_VALS <- FIT_ALL["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_ALL["xmid","Estimate"]) / FIT_ALL["scal","Estimate"] ) )
		# 	points(X_VALS,Y_VALS, type="l", col=FIT_COLS[5], lwd=6, lty=2)
		# }
		# print("Plotting Gup")
		# if (is.character(MOD_Gup.perm)==F) {
		# 	Y_VALS <- FIT_Gup.perm["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gup.perm["xmid","Estimate"]) / FIT_Gup.perm["scal","Estimate"] ) )
		# 	points(X_VALS,Y_VALS, type="l", col=FIT_COLS[4], lwd=6, lty=2)
		# }
		# print("Plotting Gdn")
		# if (is.character(MOD_Gdn.perm)==F) {
		# 	Y_VALS <- FIT_Gdn.perm["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_Gdn.perm["xmid","Estimate"]) / FIT_Gdn.perm["scal","Estimate"] ) )
		# 	points(X_VALS,Y_VALS, type="l", col=FIT_COLS[3], lwd=6, lty=2)
		# }
		# legend(.5,110, legend=c("Ind - GEx Up","Ind - GEx Down","Fit - ALL","Fit - GEx Up","Fit - GEx Down"),lty=c(1,1,2,2,2),col=FIT_COLS[c(2,1,5,4,3)], lwd=c(1,1,6,6,6) )		# if ( LL_ALL!=0 & LL.perm[p,2]!=0 & LL.perm[p,3]!=0 ) { text(-2,20, pos=4, labels=paste("p =",formatC(P_VAL_G,2,format="e")), col="black") }
		# text(.4,115, pos=4, labels="(Randomly Permuted Data)", col="black")
		# dev.off()
		COMPILE.LL[[PROBE]] <- c( LL_ALL, LL_Gup, LL_Gdn )
		COMPILE.LR[[PROBE]] <- LR_ALL_G
		COMPILE.LL.perm[[PROBE]] <- LL.perm[1:p,]
		COMPILE.LR.perm[[PROBE]] <- LR.perm[1:p]
	} # Close for loop of Probes
	COMPILE <- list( COMPILE.LL, COMPILE.LR, COMPILE.LL.perm, COMPILE.LR.perm )
	names(COMPILE) <- c("LL","LR","LL.perm","LR.perm")
	return(COMPILE)
} # Close GEXRUN function

#################################################################
## FCT: Get best candidate genes based on linear model ##########
#################################################################

## Make Function to Get Candidate Genes for a Drug
GENE_CAND <- function(drug_name, plot, pause) {
	# DOI <- "DoxorubicinHCl"
	# DOI <- "Cladribine"
	## Specify Drug Name
	DOI <- drug_name
	DOI.gex <- DRUG_KEY[ which(DRUG_KEY[,"DR"]==DOI), "GEX" ]
	print(paste("%%%%%%%% Drug =",DOI,"%%%%%%%%"))

	## Set up New Data Frame w/ Only Drug of Interest
	COLS_PULL <- c(names(DAT.2)[1:6],"Dose","iline","Batch",DOI)
	DAT.3a <- DAT.2[,COLS_PULL]
	names(DAT.3a)[ncol(DAT.3a)] <- "Resp"

	## Include all GEX data in data frame
	TEMP_GEX <- t( GEX_LIST[[DOI.gex]][,3:ncol(GEX_LIST[[DOI.gex]])] )
	TEMP_KEY <- GEX_LIST[[DOI.gex]][,1:2]
	colnames(TEMP_GEX) <- TEMP_KEY[,"Probe"]
	DAT.3b <- merge( DAT.3a, TEMP_GEX, by.x="Cell",by.y="row.names" )

	## Plot Dose-Response (linear) Fit
	if ( plot==T ) {
		par( ask=pause )
		plot( Resp ~ Dose, data=DAT.3b, pch="+", col="chartreuse2", main=PROBE.print, xlab="Dose", ylab="Response" )
		abline( lm( Resp ~ Dose, data=DAT.3b ) )
	}

	PROBE.mods <- list()
	PROBE.p <- numeric( nrow(TEMP_KEY) )
	names(PROBE.p) <- TEMP_KEY[,"Probe"]
	for ( probe_ind in 1:nrow(TEMP_KEY) ) {
		## Specify Probe to Check Out
		PROBE <- TEMP_KEY[probe_ind,"Probe"]
		PROBE.gene <- TEMP_KEY[probe_ind,"Gene"]
		PROBE.print <- paste( PROBE, " (",PROBE.gene,")",sep="" )

		## Make data frame with only 1 Probe
		DAT.3 <- data.frame(DAT.3b[,1:which(colnames(DAT.3b)=="Resp")],DAT.3b[,PROBE])
		names(DAT.3)[ncol(DAT.3)] <- "Probe"
		DAT.3 <- DAT.3[which(!is.na(DAT.3$Probe)),]
		 # Stratify Probe based on Median
		# DAT.3[,"Probe"] <- c(0,1)[ factor(DAT.3[,"Probe"] >= median(DAT.3[,"Probe"]) ) ]

		## Model Response vs GEX in linear model
		TEMP_MOD <- lm( Resp ~ Dose + Probe, data=DAT.3 )
		PROBE.mods[[PROBE]] <- TEMP_MOD
		PROBE.p[probe_ind] <- summary(TEMP_MOD)$coefficients["Probe","Pr(>|t|)"]

		## Plot GEX vs Residuals
		if ( plot==T ) {
			par( ask=pause )
			DOSE_MOD <- lm( Resp ~ Dose, data=DAT.3b )
			plot( resid(DOSE_MOD) ~ DAT.3$Probe[as.numeric(names(resid(DOSE_MOD)))], pch="+", col="deepskyblue2", main=PROBE.print, xlab="GEX", ylab="Resids" )
			if (PROBE.p[probe_ind]<.05) { COLOR <- "firebrick2" }else{ COLOR <- "black" }
			abline( lm( resid(DOSE_MOD) ~ DAT.3$Probe[as.numeric(names(resid(DOSE_MOD)))] ), col=COLOR )
		}
	}
	## Plot a few things for kicks
	EXP <- 1:length(PROBE.p) / length(PROBE.p)
	OBS <- sort( PROBE.p )
	LIM <- c(0,8)
	plot( 0,0,type="n", xlim=LIM, ylim=LIM)
	abline( 0,1 )
	points( -log10(EXP), -log10(OBS), col="deepskyblue2", pch="+" )
	## Compile Outputs
	return( PROBE.p )
} # Close GENE_CAND Function

#################################################################
## RUN IT #######################################################
############# GEXRUN(which_drug,which_list) #####################
OUT <- list()
## Run it
 # Which Drug to Run
which_drug <- "DoxorubicinHCl"
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
# OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[3:how_many_probes] )
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:2] )

## Run it
 # Which Drug to Run
which_drug <- "MLN9708MLN2238"
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:how_many_probes] )

## Run it
 # Which Drug to Run
which_drug <- "BosutinibSKI606"
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:how_many_probes] )

## Run it
 # Which Drug to Run
which_drug <- ALL_DRUGS[ grep("Criz",ALL_DRUGS) ]
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:how_many_probes] )

## Run it
 # Which Drug to Run
which_drug <- ALL_DRUGS[ grep("Palb",ALL_DRUGS) ]
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:how_many_probes] )

## Run it
 # Which Drug to Run
which_drug <- ALL_DRUGS[ grep("Dacom",ALL_DRUGS) ]
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:how_many_probes] )

## Run it
 # Which Drug to Run
which_drug <- ALL_DRUGS[ grep("Sunit",ALL_DRUGS) ]
 # Which Probes to Run
how_many_probes <- 20
GEX.p <- GENE_CAND( which_drug, plot=F, pause=F )
which_probes <- order( GEX.p )[1:how_many_probes]
OUT[[which_drug]] <- GEXRUN( which_drug, which_probes[1:how_many_probes] )






png(paste(PathToSave,"PL_1-GEX_Perm_Hist_Example.png",sep=""), width=1250,height=1050,pointsize=24)
hist( OUT$Crizotinib$LR.perm$`231666`, breaks=seq(0,160,5), col="slateblue3", main="Permuted Distribution of Likelihood Ratio Statistics", xlab="LR Stat"  )
abline( h=seq(0,15,1), lty=2, col="grey50" )
hist( OUT$Crizotinib$LR.perm$`231666`, breaks=seq(0,160,5), col="slateblue3", main="Permuted Distribution of Likelihood Ratio Statistics", xlab="LR Stat", add=T)
arrows( OUT$Crizotinib$LR$`231666`, 1, OUT$Crizotinib$LR$`231666`, 0.05, lwd=5, col="springgreen3" )
dev.off()


png(paste(PathToSave,"PL_1-GEX_Ind_Fit.png",sep=""), width=1250,height=1050,pointsize=24)
TEST <- nlsList( Resp ~ SSlogis( Dose, Asym, xmid, scal ) | Cell, data=DAT.3b, start=START, subset=which(!is.na(DAT.3b$Resp)) )
IND_COEF <- coef(TEST)
COLS.list <- c("firebrick2","chocolate2","gold2","springgreen2","steelblue2","slateblue3","black")
COLS <- colorRampPalette(COLS.list)(nrow(IND_COEF))
plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("Individual Dose-Resp - ",DOI,sep=""))
abline( h=seq(0,200,20), col="grey50", lty=c(1,2,2,2,2,1,2,2,2,2,1) )
abline( v=seq(-3,3,1), col="grey50", lty=2 )
points(DAT.3b[,"Dose"],DAT.3b[,"Resp"], pch="+", col="grey70" )
for (i in 1:nrow(IND_COEF) ) {
	Y_VALS <- IND_COEF$Asym[i] / ( 1 + exp( -(X_VALS - IND_COEF$xmid[i]) / IND_COEF$scal[i] ) )	
	points(X_VALS,Y_VALS, type="l", col=COLS[i], lwd=3 )
}
dev.off()


