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
DATE <- "20150115"

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

# #################################################################
# ## LOAD GENE EXPRESSION DATA ####################################
# #################################################################

# TOP <- GTOP <- list()

# ## Load Top Genes and Expression Data

#  # Top Gene Names
# TOP$dox <- read.table(paste(PathToTop2,"doxorubicin_top50_Genes_allcelllines.txt",sep=""),sep="\t",header=T)
# names(TOP$dox) <- c("Probe","Gene","P")
# TOP$iri <- read.table(paste(PathToTop3,"Irinotecan_top50_Genes_allcelllines.txt",sep=""),sep="\t",header=T)
# names(TOP$iri) <- c("Probe","Gene","P")
# TOP$top <- read.table(paste(PathToTop3,"Topotecan_top50_Genes_allcelllines.txt",sep=""),sep="\t",header=T)
# names(TOP$top) <- c("Probe","Gene","P")

#  # Top Expression Values
# GTOP$dox.1 <- read.table(paste(PathToTop2,"doxorubicin_top50_GEX_allcelllines.txt",sep=""),sep="\t",header=T)
# names(GTOP$dox.1)[1] <- c("Probe")
# GTOP$iri.1 <- read.table(paste(PathToTop3,"Irinotecan_top50_GEX_allcelllines.txt",sep=""),sep="\t",header=T)
# names(GTOP$iri.1)[1] <- c("Probe")
# GTOP$top.1 <- read.table(paste(PathToTop3,"Topotecan_top50_GEX_allcelllines.txt",sep=""),sep="\t",header=T)
# names(GTOP$top.1)[1] <- c("Probe")

# #################################################################
# ## ORGANIZE GENE EXPRESSION DATA ################################
# #################################################################

# ## Specify Drug Names
# ALL_DRUGS <- colnames(DAT.1)[7:134]
# BTOP <- KEY <- GEX <- KEYS <- list()
#  # and Shorthand
# DRUG_SHORT_HAND <- c("dox","iri","top")
# names(DRUG_SHORT_HAND) <- c( "DoxorubicinHCl","IrinotecanHCl","TopotecanHCl")

# ## Loop through Drugs of Interest
# for ( drug in DRUG_SHORT_HAND ) {
# 	## Create Gene-Probe Key
# 	KEY[[drug]] <- TOP[[drug]][which(!duplicated(TOP[[drug]]$Probe)),c("Gene","Probe")]
# 	KEY[[drug]] <- data.frame(KEY[[drug]],TAG=paste(KEY[[drug]]$Gene,KEY[[drug]]$Probe,sep="_"))
# 	KEY[[drug]]$TAG <- gsub("---",".",KEY[[drug]]$TAG, fixed=T)
# 	KEY[[drug]]$TAG <- gsub("-",".",KEY[[drug]]$TAG, fixed=T)
# 	KEY[[drug]]$TAG <- gsub(" ","",KEY[[drug]]$TAG, fixed=T)
# 	KEY[[drug]]$TAG <- gsub("///",".",KEY[[drug]]$TAG, fixed=T)	
	
# 	## Give Expression Values Gene Names
# 	GTOP[[drug]] <- merge(KEY[[drug]],GTOP[[paste(drug,".1",sep="")]])
# 	BTOP[[drug]] <- data.frame(GTOP[[drug]][,4:ncol(GTOP[[drug]])])
# 	for ( which_gene in 1:nrow(BTOP[[drug]]) ) {
# 		BTOP[[drug]][which_gene,] <- as.numeric( as.numeric(BTOP[[drug]][which_gene,]) > median(as.numeric(BTOP[[drug]][which_gene,])) )
# 	} # apply(GTOP[,2:ncol(GTOP)],2,median) )
# 	rownames(BTOP[[drug]]) <- GTOP[[drug]]$TAG

# 	## Compile GEX data into list
# 	GEX[[drug]] <- t( BTOP[[drug]] )
# 	KEYS[[drug]] <- data.frame(GTOP[[drug]][,1:3],PAREN=paste(GTOP[[drug]][,2]," (",GTOP[[drug]][,1],")",sep=""))
# }

#################################################################
## USE NLME TO FIT LOGISTIC MODELS ##############################
#################################################################

## Compile a few values
DOSES <- unique(DAT.2$Dose_Value)
CELL_LINES <- as.character(unique(DAT.2$Cell))

## Make Function to Compute Model & Plot this Shiz
GEXRUN <- function(drug_list) {
	for ( DOI in drug_list ) {
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

		## Run Nonlinear Fit for all cells (using xmid as Random Effect)
		print("Running MOD_ALL")
		MOD_ALL <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START ),T)
		if (is.character(MOD_ALL)==F) {
			EST_ALL <- coef(MOD_ALL)[[1]]
			SUM_ALL <- summary(MOD_ALL)
			FIT_ALL <- SUM_ALL$coefficients
			LL_ALL <- logLik(MOD_ALL)[1]
			PAR_ALL <- attr(logLik(MOD_ALL),"df")
		}else { EST_ALL <- SUM_ALL <- FIT_ALL <- list() ; LL_ALL <- PAR_ALL <- 0  }
		## Run Nonlinear Fit for BRAF Mutant Status
		 # Mutant
		print("Running MOD_BRAF1")
		MOD_BRAF1 <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=BRAF==1 ),T)
		if (is.character(MOD_BRAF1)==F) {
			EST_BRAF1 <- coef(MOD_BRAF1)[[1]]
			SUM_BRAF1 <- summary(MOD_BRAF1)
			FIT_BRAF1 <- SUM_BRAF1$coefficients
			LL_BRAF1 <- logLik(MOD_BRAF1)[1]
			PAR_BRAF1 <- attr(logLik(MOD_BRAF1),"df")
		}else { EST_BRAF1 <- SUM_BRAF1 <- FIT_BRAF1 <- list() ; LL_BRAF1 <- PAR_BRAF1 <- 0  }
		 # Wildtype
		print("Running MOD_BRAF0")
		MOD_BRAF0 <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=BRAF==-1 ),T)
		if (is.character(MOD_BRAF0)==F) {
			EST_BRAF0 <- coef(MOD_BRAF0)[[1]]
			SUM_BRAF0 <- summary(MOD_BRAF0)
			FIT_BRAF0 <- SUM_BRAF0$coefficients
			LL_BRAF0 <- logLik(MOD_BRAF0)[1]
			PAR_BRAF0 <- attr(logLik(MOD_BRAF0),"df")
		}else { EST_BRAF0 <- SUM_BRAF0 <- FIT_BRAF0 <- list() ; LL_BRAF0 <- PAR_BRAF0 <- 0  }
		## Run Nonlinear Fit for NRAS Mutant Status
		 # Mutant
		print("Running MOD_NRAS1")
		MOD_NRAS1 <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=NRAS==1 ),T)
		if (is.character(MOD_NRAS1)==F) {
			EST_NRAS1 <- coef(MOD_NRAS1)[[1]]
			SUM_NRAS1 <- summary(MOD_NRAS1)
			FIT_NRAS1 <- SUM_NRAS1$coefficients
			LL_NRAS1 <- logLik(MOD_NRAS1)[1]
			PAR_NRAS1 <- attr(logLik(MOD_NRAS1),"df")
		}else { EST_NRAS1 <- SUM_NRAS1 <- FIT_NRAS1 <- list() ; LL_NRAS1 <- PAR_NRAS1 <- 0  }
		 # Wildtype
		print("Running MOD_NRAS0")
		MOD_NRAS0 <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3, start=START, subset=NRAS==-1 ),T)
		if (is.character(MOD_NRAS0)==F) {
			EST_NRAS0 <- coef(MOD_NRAS0)[[1]]
			SUM_NRAS0 <- summary(MOD_NRAS0)
			FIT_NRAS0 <- SUM_NRAS0$coefficients
			LL_NRAS0 <- logLik(MOD_NRAS0)[1]
			PAR_NRAS0 <- attr(logLik(MOD_NRAS0),"df")
		}else { EST_NRAS0 <- SUM_NRAS0 <- FIT_NRAS0 <- list() ; LL_NRAS0 <- PAR_NRAS0 <- 0  }

		## Calculate Likelihood Ratio for 1 vs 2
		LR_ALL_BRAF <- LR_ALL_NRAS <- "NA"
		if ( LL_ALL!=0 & LL_BRAF1!=0 & LL_BRAF0!=0 ) {
			LR_ALL_BRAF <- -2*( LL_ALL - (LL_BRAF1 + LL_BRAF0) )
		}
		if ( LL_ALL!=0 & LL_NRAS1!=0 & LL_NRAS0!=0 ) {
			LR_ALL_NRAS <- -2*( LL_ALL - (LL_NRAS1 + LL_NRAS0) )
		}
		LR <- c(LR_ALL_NRAS,LR_ALL_BRAF) ; names(LR) <- c("BRAF","NRAS")

		###############################################
		## PERMUTE ## Random = ~ xmid | Cell

		## Set parameters for permutations
		n_perm <- 100
		h.lim <- 2
		LL.perm <- array( , c(n_perm,4) ) ; colnames(LL.perm) <- c("BRAF1","BRAF0","NRAS1","NRAS0")
		PAR.perm <- array( , c(n_perm,4) ) ; colnames(PAR.perm) <- c("BRAF1","BRAF0","NRAS1","NRAS0")
		LR.perm <- P.perm <- DF.perm <- array( , c(n_perm,2) ) # numeric( n_perm )
		colnames(LR.perm) <- colnames(P.perm) <- colnames(DF.perm) <- c("BRAF","NRAS")
		print("%%%%%%%% Starting Permutations %%%%%%%%")
		## Get number of cell lines w/ mutation statuses
		BRAF1.len <- length( unique( DAT.3$Cell[ which(DAT.3$BRAF==1) ] ) )
		NRAS1.len <- length( unique( DAT.3$Cell[ which(DAT.3$NRAS==1) ] ) )
		## Set up Indicator Variable for Permutations of BRAF/NRAS
		if ( LL_ALL!=0 & LL_BRAF1!=0 & LL_BRAF0!=0 ) {
			BRAF_perm <- 1
		}else{ BRAF_perm <- 0 }
		if ( LL_ALL!=0 & LL_NRAS1!=0 & LL_NRAS0!=0 ) {
			NRAS_perm <- 1
		}else{ NRAS_perm <- 0 }
		keep_permuting <- max( BRAF_perm, NRAS_perm )
		## Count Loops
		p <- p.BRAF <- p.NRAS <- 0
		start_perm <- proc.time()
		## Permute until I say stop...
		while ( keep_permuting==1 ) {
			p <- p + 1
			print(paste( "Perm:",p,"-",round(proc.time()-start_perm,3)[3] ))

			## BRAF ##
			if ( BRAF_perm==1 ) {
				print("# BRAF #")
				p.BRAF <- p.BRAF + 1
				## Create Permuted Response Table (retaining Dosing)
				DAT.3.perm <- DAT.3
				DAT.3.perm[,"BRAF"] <- -1
				WHICH_CELLS <- sample( CELL_LINES, BRAF1.len, replace=F )
				WHICH_ROWS <- which( DAT.3.perm[,"Cell"] %in% WHICH_CELLS )
				DAT.3.perm[WHICH_ROWS,"BRAF"] <- 1
				
				## Re-Run Models ##
				 # Mutant
				print("Mut")
				MOD_BRAF1.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3.perm, start=START, subset=BRAF==1 ),T)
				if (is.character(MOD_BRAF1.perm)==F) {
					EST_BRAF1.perm <- coef(MOD_BRAF1.perm)[[1]]
					SUM_BRAF1.perm <- summary(MOD_BRAF1.perm)
					FIT_BRAF1.perm <- SUM_BRAF1.perm$coefficients
					LL_BRAF1.perm <- logLik(MOD_BRAF1.perm)[1]
					PAR_BRAF1.perm <- attr(logLik(MOD_BRAF1.perm),"df")
				}else { EST_BRAF1.perm <- SUM_BRAF1.perm <- FIT_BRAF1.perm <- list() ; LL_BRAF1.perm <- PAR_BRAF1.perm <- 0  }
				 # Wildtype
				print("Wild")
				MOD_BRAF0.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3.perm, start=START, subset=BRAF==-1 ),T)
				if (is.character(MOD_BRAF0.perm)==F) {
					EST_BRAF0.perm <- coef(MOD_BRAF0.perm)[[1]]
					SUM_BRAF0.perm <- summary(MOD_BRAF0.perm)
					FIT_BRAF0.perm <- SUM_BRAF0.perm$coefficients
					LL_BRAF0.perm <- logLik(MOD_BRAF0.perm)[1]
					PAR_BRAF0.perm <- attr(logLik(MOD_BRAF0.perm),"df")
				}else { EST_BRAF0.perm <- SUM_BRAF0.perm <- FIT_BRAF0.perm <- list() ; LL_BRAF0.perm <- PAR_BRAF0.perm <- 0  }

				## Compile Permutation Results
				LL.perm[p,c("BRAF1","BRAF0")] <- c( LL_BRAF1.perm, LL_BRAF0.perm)
				PAR.perm[p,c("BRAF1","BRAF0")] <- c( PAR_BRAF1.perm, PAR_BRAF0.perm)
				LR.perm[p,"BRAF"] <- -2*( LL_ALL - (LL_BRAF1.perm + LL_BRAF0.perm) )
				 # In the LRT, degrees of freedom is equal to the number of additional parameters in the more complex model.
				DF.perm[p,"BRAF"] <- sum( PAR.perm[p,c("BRAF0","BRAF1")] )
				P.perm[p,"BRAF"] <- dchisq(LR.perm[p,"BRAF"], DF.perm[p,"BRAF"])
				## If it's clearly not going to be significant, stop permuting
				h.BRAF <- length(which( LR.perm[,"BRAF"] > LR_ALL_BRAF ))
				if ( h.BRAF >= h.lim ) { BRAF_perm <- 0 }
			}

			## NRAS ##
			if ( NRAS_perm==1 ) {
				print("# NRAS #")
				p.NRAS <- p.NRAS + 1
				## Create Permuted Response Table (retaining Dosing)
				DAT.3.perm[,"NRAS"] <- -1
				WHICH_CELLS <- sample( CELL_LINES, NRAS1.len, replace=F )
				WHICH_ROWS <- which( DAT.3.perm[,"Cell"] %in% WHICH_CELLS )
				DAT.3.perm[WHICH_ROWS,"NRAS"] <- 1
				
				## Re-Run Models ##
				 # Mutant
				print("Mut")
				MOD_NRAS1.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3.perm, start=START, subset=NRAS==1 ),T)
				if (is.character(MOD_NRAS1.perm)==F) {
					EST_NRAS1.perm <- coef(MOD_NRAS1.perm)[[1]]
					SUM_NRAS1.perm <- summary(MOD_NRAS1.perm)
					FIT_NRAS1.perm <- SUM_NRAS1.perm$coefficients
					LL_NRAS1.perm <- logLik(MOD_NRAS1.perm)[1]
					PAR_NRAS1.perm <- attr(logLik(MOD_NRAS1.perm),"df")
				}else { EST_NRAS1.perm <- SUM_NRAS1.perm <- FIT_NRAS1.perm <- list() ; LL_NRAS1.perm <- PAR_NRAS1.perm <- 0  }
				 # Wildtype
				print("Wild")
				MOD_NRAS0.perm <- try(nlmer( Resp ~ SSlogis(Dose, Asym, xmid, scal) ~ xmid | Cell, data=DAT.3.perm, start=START, subset=NRAS==-1 ),T)
				if (is.character(MOD_NRAS0.perm)==F) {
					EST_NRAS0.perm <- coef(MOD_NRAS0.perm)[[1]]
					SUM_NRAS0.perm <- summary(MOD_NRAS0.perm)
					FIT_NRAS0.perm <- SUM_NRAS0.perm$coefficients
					LL_NRAS0.perm <- logLik(MOD_NRAS0.perm)[1]
					PAR_NRAS0.perm <- attr(logLik(MOD_NRAS0.perm),"df")
				}else { EST_NRAS0.perm <- SUM_NRAS0.perm <- FIT_NRAS0.perm <- list() ; LL_NRAS0.perm <- PAR_NRAS0.perm <- 0  }

				## Compile Permutation Results
				LL.perm[p,c("NRAS1","NRAS0")] <- c( LL_NRAS1.perm, LL_NRAS0.perm)
				PAR.perm[p,c("NRAS1","NRAS0")] <- c( PAR_NRAS1.perm, PAR_NRAS0.perm)
				LR.perm[p,"NRAS"] <- -2*( LL_ALL - (LL_NRAS1.perm + LL_NRAS0.perm) )
				 # In the LRT, degrees of freedom is equal to the number of additional parameters in the more complex model.
				DF.perm[p,"NRAS"] <- sum( PAR.perm[p,] )
				P.perm[p,"NRAS"] <- dchisq(LR.perm[p], DF.perm[p])
				## If it's clearly not going to be significant, stop permuting
				h.NRAS <- length( which(LR.perm[,"NRAS"] > LR_ALL_NRAS ))
				if ( h.NRAS >= h.lim ) { NRAS_perm <- 0 }
			}
			## Check if both mutation models are done permuting
			if ( ( BRAF_perm==0 & NRAS_perm==0 ) | p == n_perm ) {
				keep_permuting <- 0
			}
		} # Close permutation "while" loop
		
		## Calculate P-Values based on Permutations
		P_PERM <- numeric(2) ; names(P_PERM) <- c("BRAF","NRAS")
		P_PERM["BRAF"] <- ( length(which( LR.perm[1:p.BRAF,"BRAF"] > LR_ALL_BRAF )) + 1 ) / (p.BRAF+1)
		P_PERM["NRAS"] <- ( length(which( LR.perm[1:p.NRAS,"NRAS"] > LR_ALL_NRAS )) + 1 ) / (p.NRAS+1)

		# par(mfrow=c(2,2))
		# plot( Resp ~ Dose, data=DAT.3, col=factor(DAT.3$BRAF) )
		# plot( Resp ~ Dose, data=DAT.3.perm, col=factor(DAT.3.perm$BRAF) )
		# plot( Resp ~ Dose, data=DAT.3, col=factor(DAT.3$NRAS) )
		# plot( Resp ~ Dose, data=DAT.3.perm, col=factor(DAT.3.perm$NRAS) )

		###########################################
		## Plot Fits based on Model Coefficients ##
		PCH_COLS.BRAF <- c("grey70","dodgerblue1")
		PCH_COLS.NRAS <- c("grey70","firebrick1")
		X_VALS <- seq(-2,2,.01)
		FIT_COLS.BRAF <- c("grey40","black","dodgerblue3")
		FIT_COLS.NRAS <- c("grey40","black","firebrick3")
		png(paste(PathToSave,"PL_1-NLME_Cell_Mut_",DOI,".png",sep=""), width=2000,height=1000,pointsize=24)
		par(mfrow=c(1,2))
		## Plot BRAF Data
		print("### Plotting BRAF Split ###")
		print("Plotting Empty Plot and MOD_ALL")
		plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("BRAF Mut & Dose-Resp - ",DOI,sep=""))
		abline( h=seq(0,200,20), col="grey50", lty=c(1,2,2,2,2,1,2,2,2,2,1) )
		abline( v=seq(-3,3,1), col="grey50", lty=2 )
		points(DAT.3[,"Dose"],DAT.3[,"Resp"], pch="+", col=PCH_COLS.BRAF[factor(DAT.3$BRAF)] )
		BRAF_COLS <- PCH_COLS.BRAF[factor(DAT.3$BRAF[seq(1,540,9)])]
		if (is.character(MOD_ALL)==F) {
			for (i in 1:nrow(EST_ALL) ) {
				Y_VALS <- EST_ALL$Asym[i] / ( 1 + exp( -(X_VALS - EST_ALL$xmid[i]) / EST_ALL$scal[i] ) )	
				points(X_VALS,Y_VALS, type="l", col=BRAF_COLS[i] )
			}
			Y_VALS <- FIT_ALL["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_ALL["xmid","Estimate"]) / FIT_ALL["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS.BRAF[2], lwd=6, lty=2)
		}
		print("Plotting Mutant")
		if (is.character(MOD_BRAF0)==F) {
			Y_VALS <- FIT_BRAF0["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_BRAF0["xmid","Estimate"]) / FIT_BRAF0["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS.BRAF[1], lwd=6, lty=2)
		}
		print("Plotting Wildtype")
		if (is.character(MOD_BRAF1)==F) {
			Y_VALS <- FIT_BRAF1["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_BRAF1["xmid","Estimate"]) / FIT_BRAF1["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS.BRAF[3], lwd=6, lty=2)
		}
		if ( LL_ALL!=0 & LL_BRAF0!=0 & LL_BRAF1!=0 ) {
			text(-2,10, pos=4, labels=paste("Permutation-Based P-Value (",p.BRAF,")",sep=""), col="black")
			text(-2,5, pos=4, labels=formatC(P_PERM["BRAF"],digits=2,format="e") ) # paste("p =",formatC(P_PERM,2,format="e")) )
		}
		## Legend
		if ( length(which(Y_VALS>60))/length(Y_VALS) > .5 ) {
			LEG.coords <- c(-2,50)
		}else{ LEG.coords <- c(.75,118) }
		legend( LEG.coords[1],LEG.coords[2], legend=c("Line - BRAF+","Line - BRAF-","Fit - ALL","Fit - BRAF+","Fit - BRAF-"),lty=c(1,1,2,2,2),col=c(PCH_COLS.BRAF,FIT_COLS.BRAF), lwd=c(1,1,6,6,6) )
		## Plot NRAS Data
		print("### Plotting NRAS Split ###")
		print("Plotting Empty Plot and MOD_ALL")
		plot(0,0,type="n", xlim=range(X_VALS), ylim=c(0,120), xlab="Concentration log10(mM)", ylab="% Cell Viability", main=paste("NRAS Mut & Dose-Resp - ",DOI,sep=""))
		abline( h=seq(0,200,20), col="grey50", lty=c(1,2,2,2,2,1,2,2,2,2,1) )
		abline( v=seq(-3,3,1), col="grey50", lty=2 )
		points(DAT.3[,"Dose"],DAT.3[,"Resp"], pch="+", col=PCH_COLS.NRAS[factor(DAT.3$NRAS)] )
		NRAS_COLS <- PCH_COLS.NRAS[factor(DAT.3$NRAS[seq(1,540,9)])]
		if (is.character(MOD_ALL)==F) {
			for (i in 1:nrow(EST_ALL) ) {
				Y_VALS <- EST_ALL$Asym[i] / ( 1 + exp( -(X_VALS - EST_ALL$xmid[i]) / EST_ALL$scal[i] ) )	
				points(X_VALS,Y_VALS, type="l", col=NRAS_COLS[i] )
			}
			Y_VALS <- FIT_ALL["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_ALL["xmid","Estimate"]) / FIT_ALL["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS.NRAS[2], lwd=6, lty=2)
		}
		print("Plotting Mutant")
		if (is.character(MOD_NRAS0)==F) {
			Y_VALS <- FIT_NRAS0["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_NRAS0["xmid","Estimate"]) / FIT_NRAS0["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS.NRAS[1], lwd=6, lty=2)
		}
		print("Plotting Wildtype")
		if (is.character(MOD_NRAS1)==F) {
			Y_VALS <- FIT_NRAS1["Asym","Estimate"] / ( 1 + exp( -(X_VALS - FIT_NRAS1["xmid","Estimate"]) / FIT_NRAS1["scal","Estimate"] ) )
			points(X_VALS,Y_VALS, type="l", col=FIT_COLS.NRAS[3], lwd=6, lty=2)
		}
		if ( LL_ALL!=0 & LL_NRAS0!=0 & LL_NRAS1!=0 ) {
			text(-2,10, pos=4, labels=paste("Permutation-Based P-Value (",p.NRAS,")",sep=""), col="black")
			text(-2,5, pos=4, labels=formatC(P_PERM["NRAS"],digits=2,format="e") ) # paste("p =",formatC(P_PERM,2,format="e")) )
		}
		## Legend
		legend( LEG.coords[1],LEG.coords[2], legend=c("Line - NRAS+","Line - NRAS-","Fit - ALL","Fit - NRAS+","Fit - NRAS-"),lty=c(1,1,2,2,2),col=c(PCH_COLS.NRAS,FIT_COLS.NRAS), lwd=c(1,1,6,6,6) )
		dev.off()
	}
} # Close GEXRUN function

#################################################################
## RUN IT #######################################################
############# GEXRUN(drug_name,start,stop) ######################

## Potential Drugs of Interest
# DRUGS.pref <- c("Trametinib","MEK162","Palbociclib","MLN0128","GSK2141795")
DRUG_LIST <- names(DAT.2)[7:134]
DRUGS.pref2 <- c("Trametinib","MEK162","PD325901","Palbociclib","MLN0128","INK128","OSI027","Everolimus","Sirilimus","Temsirilimus","GSK2141795")
DRUGS.ccle <- c("Topotecan","Nilotinib","Lapatinib","Irinotecan","Erlotinib","Sorafenib")
DRUGS.prev <- DRUG_LIST[ c(4,11,17,18,19,21,23,33,41,47,49,64,66,69,73,79,96,100,101,112,119,122,125,126) ]
DRUGS.fit <- DRUG_LIST[ c(11,18,33,47,69,122,126) ]
DRUGS.gen <- c("Cabozantinib","Dacomitinib","Etoposide","MLN2480","Palbociclib","Sunitinib","Vorinostat")

## Get actual names of Drugs from DRUGS.pref
# drug_list <- c()
# for (drug in DRUGS.pref) {
# 	drug_list <- c( drug_list, DRUG_LIST[grep(drug,DRUG_LIST)] )
# }
# GEXRUN(drug_list)

## Get actual names of Drugs from DRUGS.ccle
drug_list <- c()
for (drug in DRUGS.ccle) {
	drug_list <- c( drug_list, DRUG_LIST[grep(drug,DRUG_LIST)] )
}
GEXRUN(drug_list)

## Get actual names of Drugs from DRUGS.pref2
drug_list <- c()
for (drug in DRUGS.pref2) {
	drug_list <- c( drug_list, DRUG_LIST[grep(drug,DRUG_LIST)] )
}
GEXRUN(drug_list)

## Run for list of good fits
GEXRUN(DRUGS.fit)

## Run for previous list (i think of good fits)
DONE <- Reduce( union, list(DRUGS.fit, DRUGS.pref2, DRUGS.ccle) )
drug_list <- setdiff( DRUGS.prev, DONE )
GEXRUN(drug_list)

## Get actual names of Drugs from DRUGS.gen
drug_list <- c()
for (drug in DRUGS.gen) {
	drug_list <- c( drug_list, DRUG_LIST[grep(drug,DRUG_LIST)] )
}
GEXRUN(drug_list)


# GEXRUN("DoxorubicinHCl")
# GEXRUN("Dacarbazine")
# GEXRUN("Sorafenib")
# GEXRUN("Etoposide")
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



