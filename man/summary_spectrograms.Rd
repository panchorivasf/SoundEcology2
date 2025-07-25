% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_spectrograms.R
\name{summary_spectrograms}
\alias{summary_spectrograms}
\title{Generate and Save Summary Spectrograms}
\usage{
summary_spectrograms(
  summary_df,
  parent_dir = NULL,
  output_dir = "summary_spectros_se2",
  rmdcoff = TRUE,
  norm = FALSE
)
}
\arguments{
\item{summary_df}{A data frame containing summary statistics of recordings. It should include
columns \code{sensor_id}, \code{channel}, \code{index}, and columns named \verb{closest_to_<stat>} for categories \code{mean},
\code{max}, \code{min}, \code{median}, and \code{mode}, each referring to the closest wave file for each statistic.}

\item{parent_dir}{A character string specifying the parent directory containing the \code{.wav} files.
Defaults to \code{NULL}, which sets \code{parent_dir} to the current working directory.}

\item{output_dir}{A character string specifying the directory to save generated spectrogram images.
aults to \code{"summary_spectrograms_se2"}.}

\item{rmdcoff}{A logical value indicating whether to remove the DC offset from wave files before
generating spectrograms. Defaults to \code{TRUE}.}

\item{norm}{A logical value indicating whether to normalize the spectrogram. Defaults to \code{FALSE}.}
}
\value{
None. Spectrogram images are saved in the specified output directory.
}
\description{
This function generates and saves spectrogram images for various summary
statistics (mean, max, min, median, mode) for each recording in the input
data frame. Spectrograms are created for specific wave files as indicated in
the data, with options to normalize and remove DC offset.
}
\examples{
\dontrun{
Example usage with a data frame `summary_df` that includes necessary summary columns 
summary_spectrograms(summary_df, parent_dir = "path/to/wav/files", output_dir = "spectrograms", 
rmdcoff = TRUE, norm = FALSE)
}
}
