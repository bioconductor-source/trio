\name{combineNtrios}
\alias{combineNtrios}
\alias{combineNtriosList}


\title{
Combining Matrices with Frequencies of Genotype Combinations
}
\description{
Combines, i.e. adds together, two or a list of matrices generated with \code{colNtrios} that contain the frequencies of the ten different genotype combinations a SNP can show for an additive effect in a case-parent trio study.
}
\usage{
combineNtrios(matNtrios1, matNtrios2, idSNP1 = NULL, idSNP2 = NULL, showInfo = TRUE, 
   orderSNPs = TRUE, idAsRownames = FALSE)
	
combineNtriosList(list.Ntrios, list.idSNP = NULL, showInfo = TRUE, orderSNPs = TRUE, 
   idAsRownames = FALSE)
}

\arguments{
  \item{matNtrios1}{a numeric matrix or list of two matrices containing the numbers of trios showing the different combinations of genotypes of the offsprings and their parents, i.e. the output of \code{\link{colNtrios}}, where the arguments \code{ref} and \code{alt} must have been specified when generating \code{matNtrios1} with \code{\link{colNtrios}}. 
}
  \item{matNtrios2}{a second numeric matrix or list of two matrices that should be combined with \code{matNtrios1}, where each matrix contains the numbers of trios showing the different combinations of genotypes of the offsprings and their parents, i.e. the output of \code{\link{colNtrios}}. In the generation of \code{matNtrios2} with \code{\link{colNtrios}}, the arguments \code{ref} and \code{alt} must have been specified.  
}
  \item{idSNP1}{a character vector containing IDs for the SNPs, e.g., their names, in \code{matNtrios1}. If a SNP is available in both \code{matNtrios1} and \code{matNtrios2}, then the ID for this SNP in \code{idSNP1} and \code{idSNP2} must be exactly the same. Otherwise, this SNP is treated as two different SNPs when combining \code{matNtrios1} and \code{matNtrios2}. If not specified, the row names of \code{matNtrios1} will be used as IDs for the SNPs.
}
  \item{idSNP2}{a character vector containing IDs for the SNPs, e.g., their names, in \code{matNtrios2}. If a SNP is available in both \code{matNtrios1} and \code{matNtrios2}, then the ID for this SNP in \code{idSNP1} and \code{idSNP2} must be exactly the same. Otherwise, this SNP is treated as two different SNPs when combining \code{matNtrios1} and \code{matNtrios2}. If not specified, the row names of \code{matNtrios2} will be used as IDs for the SNPs.
}
  \item{showInfo}{logical specifying whether information on how many SNPs appear only in one of the matrices and how many appear in both matrices.
}
  \item{orderSNPs}{logical whether the SNPs should be ordered alphanumerically (by their row name) in the combined matrix. If \code{FALSE}, the rows of combined matrix contain first the SNPs that are only in \code{matNtrio1}, then the SNPs that are in both matrices, and last the SNPs that are only in \code{matNtrio1}. 
}
  \item{idAsRownames}{logical speciying whether the SNP IDs specified by \code{idSNP1} and \code{idSNP2} are used as row/SNP names in the combined matrix. If \code{FALSE} or if \code{idSNP1} and \code{idSNP2} are not specified, the row/SNP names of \code{matNtrios1} and \code{matTrios2} are also used in the combined matrix. 
}

  \item{list.Ntrios}{a list of numeric matrices or lists of two matrices that should be combined, where each matrix contains the numbers of trios showing the different combinations of genotypes of the offsprings and their parents. Thus, each element of \code{list.Ntrios} must be an output of \code{\link{colNtrios}}, where in the generation of these elements the arguments \code{ref} and \code{alt} of \code{\link{colNtrios}} must have been specified. 
}
	
  \item{list.idSNP}{a list of the same length as \code{list.Ntrios} consisting of character vector that contain the IDs for the SNPs, e.g., their names, in the corresponding element/matrix in \code{list.Ntrios}. If a SNP is available in two or more of the matrices in \code{list.Ntrios}, then the ID for this SNP must be exactly the same in all corresponding elements of \code{list.idSNP}. Otherwise, this SNP is treated as different SNPs when combining the matrices. If not specified, the row names of the matrices in \code{list.idSNP} will be used as IDs for the SNPs.
}
}

\details{
\code{combineNtriosList} is a wrapper for applying \code{combineNtrios} \eqn{n - 1} times combining iteratively the \eqn{n} elements of \code{list.Ntrios}.
  
If a SNP is available in both \code{mat.Ntrio1} and \code{mat.Ntrio2} or in two or more of the matrices in \code{list.Ntrios}, then the ID for this SNP must be exactly the same in all matrices to ensure that this SNP is identified as the same SNP in all matrices. If this is not the case, the SNP is treated as if it were a different SNP in each matrix in which the ID differs. 

If the reference allele and alternative allele of a SNP are interchanged between matrices, then \code{combineNtrios} will take this into account by exchanging the reference and the alternative allele of this SNP in \code{matNtrios2} (or the later element in \code{list.Ntrios}) and by swapping the numbers for the different genotype combinations in \code{matNtrios2} correspondingly.
} 

\value{
Depending on the specification of the argument \code{env} in the calls of \code{\link{colNtrios}} with which the input matrices of \code{combineNtrios} or \code{combineNtriosList} have been generated, a matrix or a list of two matrices containing the combination of the input matrices, where the entries for the SNPs available in two or more input matrices are combined by adding the numbers in these input matrices.

For details on the input matrices/lists, i.e. the output of \code{\link{colNtrios}}, and therefore, on the output matrix/list of \code{combineNtrios} or \code{combineNtriosList}, see \code{\link{colNtrios}}.
}

\references{
Schwender, H., Taub, M.A., Beaty, T.H., Marazita, M.L., and Ruczinski, I. (2011). 
Rapid Testing of SNPs and Gene-Environment Interactions in Case-Parent Trio Data Based on 
Exact Analytic Parameter Estimation. \emph{Biometrics}, 68, 766-773. 
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{colNtrios}}
}

\examples{
# Load the simulated trio data.
data(trio.data)

# Specify reference and alternative alleles
# of the six SNPs in mat.test, where the
# two alleles of a SNP must differ.
ref <- rep(c("A", "T", "G"), 2)
alt <- c("T", "G", "A", "C", "C", "C")

# Let us assume that the following two matrices
# contain trio data from two different
# case-parent trio studies.
mat.trio1 <- mat.test[ , 1:4]
mat.trio2 <- mat.test[ , 3:6]

# Compute the numbers of trios showing
# the different genotype combinations
# for both matrices with trio data.
mat.ntrios1 <- colNtrios(mat.trio1, ref = ref[1:4], 
   alt = alt[1:4])
mat.ntrios2 <- colNtrios(mat.trio2, ref = ref[3:6], 
   alt = alt[3:6])

# Then, these matrices can be combined by
combineNtrios(mat.ntrios1, mat.ntrios2)

# For comparison
colNtrios(mat.test)  
}

\keyword{array}
\keyword{manip}
