\name{ChemmineDrugs}
\title{
	Chemmine Drugs
}
\description{
	This package contains pre-built SQLite database of the DUD and DrugBank data sets.
	These databases can be used with functions in the ChemmineR package to perform compound 
	searches and also to retrieve pre-computed features for each compound.
}
\author{
	Kevin Horan
}
\seealso{
	\code{\link{ChemmmineR}}
}

\examples{
	library(ChemmineDrugs)
	## list included databases
	dir(system.file("extdata",package="ChemmineDrugs"))
}
