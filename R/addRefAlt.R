addRefAlt <- function(matNumber, ref, alt){
	checkMatN(matNumber)
	n.snp <- if(is.list(matNumber)) nrow(matNumber$ntriosEnv0) 
		else nrow(matNumber)
	checkRefAlt(ref, alt, n.snp)
	if(!is.null(attributes(matNumber)$REF))
		stop("REF and ALT are already available in matNumber.")
	attr(matNumber, "REF") <- ref
	attr(matNumber, "ALT") <- alt
	if(is.list(matNumber))
		class(matNumber) <- "tabNtriosGxE"
	else
		class(matNumber) <- "tabNtrios"
	matNumber
}


checkRefAlt <- function(ref, alt, nsnp){
	if(is.null(ref) | is.null(alt))
		stop("Both ref and alt must be specified.")
	if(length(ref) != nsnp)
		stop("The length of ref must be equal to the number of SNPs.")
	if(length(alt) != nsnp)
		stop("The length of alt must be equal to the number of SNPs.")
	if(any(is.na(ref)))
		stop("No missing values allowed in ref.")
	if(any(is.na(alt)))
		stop("No missing values allowed in alt.")
	if(!all(ref %in% c("A", "T", "C", "G")))
		stop("At least one value in ref is not A, T, C, or G.")
	if(!all(alt %in% c("A", "T", "C", "G")))
		stop("At least one value in alt is not A, T, C, or G.")
	if(any(ref == alt))
		stop("At least one entry/allel in ref is identical to the corresponding\n",
		 "  entry/allel in alt.")
}
