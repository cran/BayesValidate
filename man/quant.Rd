\name{quant}
\alias{quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate empirical quantile of the first entry in a vector }
\description{
  \code{quant} inputs a vector and returns the empirical quantile of the first 
	argument in the vector with respect to all entries in the vector.  Used as 
	part of the function validate for Bayesian software validation, this function 
	is used to calculate the empirical quantile of a "true" parameter value with 
	respect to a collection of posterior draws of that parameter.
}
\usage{
quant(draws)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{draws}{ Vector of parameter draws, with entry of interest, i.e., the 
	value whose quantile is being calculated, at the beginning. }
}
\details{
  Calculates the rank of the first entry of the vector with respect to all other entries, 
	subtracts .5, and divides by the length of the vector.  
}
\value{
  The empirical quantile of the first entry of the vector, a
	scalar between 0 and 1.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Samantha Cook \email{cook@stat.columbia.edu} }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
	set.seed(314)
	x<-rnorm(1000)
	quant(x)

}
\keyword{distribution}