combineNtrios <- function(matNtrios1, matNtrios2, idSNP1 = NULL, idSNP2 = NULL, showInfo = TRUE,
		orderSNPs = TRUE, idAsRownames = FALSE){
	isGxE <- checkMatNumbers(matNtrios1, matNtrios2)
	if(is.null(attributes(matNtrios1)$REF))
		stop("No information on the reference and alternative alleles of the SNPs available\n",
			"  in matNtrios1.")
	if(is.null(attributes(matNtrios2)$REF))
		stop("No information on the reference and alternative alleles of the SNPs available\n",
			"  in matNtrios2.")
	if(is.null(idSNP1) & is.null(idSNP2)){
		cat("NOTE: Since idSNP1 and idSNP2 have not been specified, the rownames of\n",
			"matNtrios1 and matNtrios2 are used for matching the SNPs.\n\n", sep = "")
		idSNP1 <- if(isGxE) rownames(matNtrios1$ntriosEnv0) else rownames(matNtrios1)
		idSNP2 <- if(isGxE) rownames(matNtrios2$ntriosEnv0) else rownames(matNtrios2)
	}
	else{
		if(is.null(idSNP1) | is.null(idSNP2))
			stop("Either both or none of idSNP1 and idSNP2 should be specified.")
		nsnps1 <- if(isGxE) nrow(matNtrios1$ntriosEnv0) else nrow(matNtrios1)
		nsnps2 <- if(isGxE) nrow(matNtrios2$ntriosEnv0) else nrow(matNtrios2)
		if(length(idSNP1) != nsnps1)
			stop("The length of idSNP1 differs from the number of SNPs in matNtrios1.")
		if(length(idSNP2) != nsnps2)
			stop("The length of idSNP2 differs from the number of SNPs in matNtrios2.")	
	}	
	if(isGxE)
		out <- combineTabNtriosGxE(matNtrios1, matNtrios2, idSNP1, idSNP2, showInfo = showInfo,
			orderSNPs = orderSNPs, idAsRownames = idAsRownames)
	else
		out <- combineTabNtrios(matNtrios1, matNtrios2, idSNP1, idSNP2, showInfo = showInfo,
			orderSNPs = orderSNPs, idAsRownames = idAsRownames)
	out
}

combineNtriosList <- function(list.Ntrios, list.idSNP = NULL, showInfo = TRUE, orderSNPs = TRUE, 
		idAsRownames = FALSE){
	if(!is.list(list.Ntrios))
		stop("list.Ntrios must be a list.")
	n.mat <- length(list.Ntrios)
	if(n.mat < 2)
		stop("list.Ntrios must consist of at least two elements.")
	if(!is.null(list.idSNP)){
		if(!is.list(list.idSNP))
			stop("list.idSNP must be a list.") 
		if(length(list.idSNP) != n.mat)
			stop("The lengths of list.idSNP and list.Ntrios differ.")
	}
	else
		list.idSNP <- vector("list", n.mat)
	matSum <- list.Ntrios[[1]]
	for(i in 2:n.mat){
		if(showInfo){
			if(i == 2)
				cat("Combination of first two matrices:\n\n", sep = "")
			else
				cat("Combination with ", i, "-th matrix:\n\n", sep = "")
		}	
		matSum <- combineNtrios(matSum, list.Ntrios[[i]], idSNP1 = list.idSNP[[i-1]],
			idSNP2 = list.idSNP[[i]], showInfo = showInfo, orderSNPs = orderSNPs,
			idAsRownames = idAsRownames)
		if(showInfo && i != n.mat)
				cat("\n", sep = "")
	}
	matSum
}


combineTabNtrios <- function(matN1, matN2, idSNP1, idSNP2, showInfo = TRUE, orderSNPs = TRUE,
		idAsRownames = FALSE){
	isUni1 <- !idSNP1 %in% idSNP2
	isUni2 <- !idSNP2 %in% idSNP1
	n.uni1 <- sum(isUni1)
	n.uni2 <- sum(isUni2)
	n.both <- sum(!isUni1)
	n.snps1 <- length(idSNP1)
	n.snps <- n.snps1 + n.uni2
	matComb <- matrix(nrow = n.snps, ncol = ncol(matN1))
	if(n.uni1 > 0)
		matComb[1:n.uni1, ] <- matN1[isUni1, ]
	if(n.uni2 > 0)
		matComb[(n.snps1 + 1): n.snps, ] <- matN2[isUni2, ]
	if(n.both > 0)
		matComb[(n.uni1 + 1):n.snps1, ] <- addTabNtrios(matN1, matN2, idSNP1, idSNP2, 
			isUni1, isUni2)
	if(showInfo)
		cat("SNPs only in first matrix: ", n.uni1, "\n", "SNPs only in second matrix: ", n.uni2,
			"\n", "SNPs in both matrices: ", n.both, "\n\n", sep = "")
	if(idAsRownames)
		rn <- c(idSNP1[isUni1], idSNP1[!isUni1], idSNP2[isUni2])
	else{
		rn1 <- rownames(matN1)
		rn <- c(rn1[isUni1], rn1[!isUni1], rownames(matN2)[isUni2])
	}
	rownames(matComb) <- rn
	colnames(matComb) <- colnames(matN1)
	ref1 <- attributes(matN1)$REF
	ref <- c(ref1[isUni1], ref1[!isUni1], attributes(matN2)$REF[isUni2]) 
	alt1 <- attributes(matN1)$ALT
	alt <- c(alt1[isUni1], alt1[!isUni1], attributes(matN2)$ALT[isUni2]) 
	if(orderSNPs){
		orn <- order(rn)
		matComb <- matComb[orn, ]
		ref <- ref[orn]
		alt <- alt[orn]
	}
	attr(matComb, "REF") <- ref
	attr(matComb, "ALT") <- alt
	class(matComb) <- "tabNtrios"
	matComb
}

addTabNtrios <- function(matN1, matN2, idSNP1, idSNP2, isUni1, isUni2){
	ref1 <- attributes(matN1)$REF[!isUni1]
	ref2 <- attributes(matN2)$REF[!isUni2]
	alt1 <- attributes(matN1)$ALT[!isUni1]
	alt2 <- attributes(matN2)$ALT[!isUni2]
	matN1 <- matN1[!isUni1, ]
	matN2 <- matN2[!isUni2, ]
	idSNP1 <- idSNP1[!isUni1]
	idSNP2 <- idSNP2[!isUni2]
	idM <- match(idSNP2, idSNP1)
	matN2 <- matN2[idM, ]
	ref2 <- ref2[idM]
	alt2 <- alt2[idM]
	sameRef <- ref1 == ref2
	if(any(alt1[sameRef] != alt2[sameRef]))
		stop("For at least one SNP present in both matrices that has the same reference\n",
			"  allele in both matrices, the alternative alleles differ.")
	if(any(ref1[!sameRef] != alt2[!sameRef]) | any(ref2[!sameRef] != alt1[!sameRef]))
		stop("For at least one SNP present in both matrices, the two alleles differ that\n",
			"  this SNP show differ between the two matrices.")
	if(any(!sameRef))
		matN2[!sameRef, ] <- revMatNumber(matN2[!sameRef, ])
	 matN1 + matN2
}

combineTabNtriosGxE <- function(matN1, matN2, idSNP1, idSNP2, showInfo = TRUE, orderSNPs = TRUE,
		idAsRownames = FALSE){
	isUni1 <- !idSNP1 %in% idSNP2
	isUni2 <- !idSNP2 %in% idSNP1
	n.uni1 <- sum(isUni1)
	n.uni2 <- sum(isUni2)
	n.both <- sum(!isUni1)
	n.snps1 <- length(idSNP1)
	n.snps <- n.snps1 + n.uni2
	matComb0 <- matComb1 <- matrix(nrow = n.snps, ncol = ncol(matN1$ntriosEnv0))
	if(n.uni1 > 0){
		matComb0[1:n.uni1, ] <- matN1$ntriosEnv0[isUni1, ]
		matComb1[1:n.uni1, ] <- matN1$ntriosEnv1[isUni1, ]
	}
	if(n.uni2 > 0){
		matComb0[(n.snps1 + 1): n.snps, ] <- matN2$ntriosEnv0[isUni2, ]
		matComb1[(n.snps1 + 1): n.snps, ] <- matN2$ntriosEnv1[isUni2, ]
	}
	if(n.both > 0){
		tmp <- addTabNtriosGxE(matN1, matN2, idSNP1, idSNP2, isUni1, isUni2)
		matComb0[(n.uni1 + 1):n.snps1, ] <- tmp$mat0
		matComb1[(n.uni1 + 1):n.snps1, ] <- tmp$mat1
	}
 
	if(showInfo)
		cat("SNPs only in first matrix list: ", n.uni1, "\n", "SNPs only in second matrix list: ", n.uni2,
			"\n", "SNPs in both matrix lists: ", n.both, "\n\n", sep = "")
	if(idAsRownames)
		rn <- c(idSNP1[isUni1], idSNP1[!isUni1], idSNP2[isUni2])
	else{
		rn1 <- rownames(matN1$ntriosEnv0)
		rn <- c(rn1[isUni1], rn1[!isUni1], rownames(matN2$ntriosEnv0)[isUni2])
	}
	rownames(matComb0) <- rownames(matComb1) <- rn
	colnames(matComb0) <- colnames(matComb1) <- colnames(matN1$ntriosEnv0)
	ref1 <- attributes(matN1)$REF
	ref <- c(ref1[isUni1], ref1[!isUni1], attributes(matN2)$REF[isUni2]) 
	alt1 <- attributes(matN1)$ALT
	alt <- c(alt1[isUni1], alt1[!isUni1], attributes(matN2)$ALT[isUni2]) 
	if(orderSNPs){
		orn <- order(rn)
		matComb0 <- matComb0[orn, ]
		matComb1 <- matComb1[orn, ]
		ref <- ref[orn]
		alt <- alt[orn]
	}
	listComb <- list(ntriosEnv0 = matComb0, ntriosEnv1 = matComb1)
	attr(listComb, "REF") <- ref
	attr(listComb, "ALT") <- alt
	class(listComb) <- "tabNtriosGxE"
	listComb
}

addTabNtriosGxE <- function(matN1, matN2, idSNP1, idSNP2, isUni1, isUni2){
	ref1 <- attributes(matN1)$REF[!isUni1]
	ref2 <- attributes(matN2)$REF[!isUni2]
	alt1 <- attributes(matN1)$ALT[!isUni1]
	alt2 <- attributes(matN2)$ALT[!isUni2]
	mat0N1 <- matN1$ntriosEnv0[!isUni1, ]
	mat1N1 <- matN1$ntriosEnv1[!isUni1, ]
	mat0N2 <- matN2$ntriosEnv0[!isUni2, ]
	mat1N2 <- matN2$ntriosEnv1[!isUni2, ]
	idSNP1 <- idSNP1[!isUni1]
	idSNP2 <- idSNP2[!isUni2]
	idM <- match(idSNP2, idSNP1)
	mat0N2 <- mat0N2[idM, ]
	mat1N2 <- mat1N2[idM, ]
	ref2 <- ref2[idM]
	alt2 <- alt2[idM]
	sameRef <- ref1 == ref2
	if(any(alt1[sameRef] != alt2[sameRef]))
		stop("For at least one SNP present in both matrices that has the same reference\n",
			"  allele in both matrices, the alternative alleles differ.")
	if(any(ref1[!sameRef] != alt2[!sameRef]) | any(ref2[!sameRef] != alt1[!sameRef]))
		stop("For at least one SNP present in both matrices, the two alleles differ that\n",
			"  this SNP show differ between the two matrices.")
	if(any(!sameRef)){
		mat0N2[!sameRef, ] <- revMatNumber(mat0N2[!sameRef, ])
		mat1N2[!sameRef, ] <- revMatNumber(mat1N2[!sameRef, ])
	}
	list(mat0 = mat0N1 + mat0N2, mat1 = mat1N1 + mat1N2)
}
	
	
checkMatNumbers <- function(matN1, matN2){
	type1 <- class(matN1)[1]
	type2 <- class(matN2)[1]
	if(!type1 %in% c("tabNtrios", "tabNtriosGxE"))
		stop("matNtrios1 does not seem to be the output of colNtrios.")
	if(!type2 %in% c("tabNtrios", "tabNtriosGxE"))
		stop("matNtrios2 does not seem to be the output of colNtrios.")
	if(type1 != type2)
		stop("One of matNtrios1 and matNtrios2 seems to be generated with the argument env\n",
			"  of colNtrios specified, while the other was not.")
	isGxE <- type1 == "tabNtriosGxE"
	if(isGxE){
		if(ncol(matN1$ntriosEnv0) != ncol(matN2$ntriosEnv1))
			stop("The number of columns differs between matNtrios1 and matNtrios2.")
	}
	else{
		if(ncol(matN1) != ncol(matN2))
			stop("The number of columns differs between matNtrios1 and matNtrios2.")
	}
	isGxE
}		
	

revMatNumber <- function(matNumber){
	if(is.list(matNumber)){
		out <- lapply(matNumber, revMatNumber)
		attr(out, "REF") <- attributes(matNumber)$ALT
		attr(out, "ALT") <- attributes(matNumber)$REF
		class(out) <- "tabNtriosGxE"
		return(out)
	}
	idRev <- c(4, 3, 2, 1, 7, 6, 5, 8, 10, 9)[1:ncol(matNumber)]
	matRev <- matNumber[ , idRev]
	colnames(matRev) <- colnames(matNumber)
	attr(matRev, "REF") <- attributes(matNumber)$ALT
	attr(matRev, "ALT") <- attributes(matNumber)$REF
	class(matRev) <- "tabNtrios"
	matRev
} 
