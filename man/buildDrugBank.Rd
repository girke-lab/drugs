\name{buildDrugBank}
\alias{buildDrugBank}
\title{
Build Drug Bank
}
\description{
	Download the Drug Bank data set and import it into an
	SQLite database. Several features are also computed, stored,
	and indexed in the database.

	This function is not needed to make use of the data sets contained in this
	package. It can be used to re-generate the dataset if desired however.
}
\usage{
buildDrugBank(dbName, downloadDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbName}{
		The name of the SQLite database to produce.
	}
  \item{downloadDir}{
		The directory into which the Drug Bank data files will be downloaded. This directory
		could be removed after the database is built.
	}
}
\value{
	The compound id numbers of the loaded compounds are returned as a vector.
}
\author{
	Kevin Horan
}

\examples{
	\dontrun{
		loadedIds = buildDrugBank("drugbank.db","dataSrc/drugbank")
	}
}

