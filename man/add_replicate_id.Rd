% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_replicate_id.R
\name{add_replicate_id}
\alias{add_replicate_id}
\title{Add replicate id for a dataset with severall treatments}
\usage{
add_replicate_id(data, treatment = "treatment")
}
\arguments{
\item{data}{a data frame}

\item{treatment}{the name of the column containing the treatment identifiers}
}
\value{
a data frame with treatment_id column
}
\description{
Add replicate id for a dataset with severall treatments
}
\examples{
add_replicate_id(my_dataset, treatment = 'zone')
}
