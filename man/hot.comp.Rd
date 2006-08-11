\name{hot.comp}
\alias{hot.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hotelling T2 }
\description{
Computes Hotelling \eqn{T^2} statistics
}
\usage{
hot.comp(expr.dat, status, set.size, nset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr.dat}{ Expression data: a matrix with genes in rows and
    individuals in columns}
  \item{status}{Label of individuals(0 or 1) with length the number of
    column in \code{expr.dat}}
  \item{set.size}{Number of genes used to discriminate between groups:
    an integer (between 2 and 6). }
  \item{nset}{Number of set of genes for which the \eqn{T^2} statistic is stored in the
    returned matrix: an integer (should be \eqn{> 1})}
}
\details{}

\value{A matrix with \code{set.size}+2 columns and \code{nset} lines.
  The index of genes are stored in columns 1 to \code{set.size}.
  Column \code{set.size}+1 contains the \eqn{T^2} statistics.
Column \code{set.size}+2 contains the corresponding nominal p-value.}

\references{Discrimination and scoring  using small sets of genes for 
    two-sample microarray data,
    G. Guillot, M. Olsson, M. Bensson, M. Rudemo. To appear in
    Mathematical Biosciences }

\author{ Gilles Guillot}
\note{}
\examples{

\dontrun{
## Assuming expression data are stored in a matrix named expr.dat
## with 10 controls + 10 patients
## the following instruction
res <- hot.comp(expr.dat=expr.dat,
                status=c(rep(0,10),rep(1,10)),
                set.size=4,
                nset=100)
## will return the 100 sets of 4 genes with highest T^2 statistic
}

}
\keyword{}
%\keyword{Expression data, Hotelling \eqn{T^2}}

