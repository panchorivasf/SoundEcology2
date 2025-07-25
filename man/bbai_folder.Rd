% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbai_folder.R
\name{bbai_folder}
\alias{bbai_folder}
\title{Calculate BBAI for all the Files in a Folder}
\usage{
bbai_folder(
  folder = NULL,
  recursive = FALSE,
  list = NULL,
  start = 0,
  end = 1,
  unit = "minutes",
  channel = "each",
  save.csv = TRUE,
  csv.name = "bbai_results.csv",
  hpf = 0,
  freq.res = 50,
  cutoff = -60,
  click.length = 10,
  difference = 10,
  gap.allowance = 2,
  spectrogram = FALSE,
  dark.plot = FALSE,
  plot.title = "",
  verbose = FALSE,
  n.cores = -1
)
}
\arguments{
\item{folder}{Character. The path to a folder with the wave files to analyze.}

\item{recursive}{Logical. Whether to search in subfolders. Default is TRUE.}

\item{list}{An optional list (subset) of files in the folder to analyze. If provided,
files outside the list will be excluded.}

\item{start}{numerical. Where to start reading the Wave.}

\item{end}{numerical. Where to end reading the Wave.}

\item{unit}{character. Unit of measurement for 'start' and 'end'. Options are
'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.}

\item{channel}{Character. If Wave is stereo and you want to use only one channel, pass either
"left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix".
If NULL (default), results are returned for each channel.}

\item{save.csv}{logical. Whether to save a csv in the working directory.}

\item{csv.name}{Character. Name for the csv file. Default is "bbai_results.csv".}

\item{hpf}{Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless
signals of interest are below that threshold.}

\item{freq.res}{Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency
bin and, therefore, the window length to be used (sampling rate / frequency resolution).}

\item{cutoff}{Numeric. The amplitude threshold (in dBFS) for removing low-amplitude values in the spectrogram.
Default is \code{-50}.}

\item{difference}{Numeric. The maximum difference in amplitude between adjacent frequency bins to be
considered part of a single 'click'. Default is \code{20}.}

\item{gap.allowance}{Numeric. The size of gaps (in frequency bins) allowed between contiguous parts of a click.
Default is \code{2}. Gaps larger than this value will split clicks.}

\item{spectrogram}{Logical. Should a spectrogram with highlighted clicks be plotted? Default is \code{TRUE}.}

\item{dark.plot}{Logical. Should the plot use a dark theme (black background)? Default is \code{FALSE}.}

\item{plot.title}{Character. The title for the plot, if \code{plot} is \code{TRUE}. Default is \code{NULL}.}

\item{verbose}{Logical. If TRUE, details of dynamic range will be printed on the console.}

\item{n.cores}{Numeric. Number of cores to be used in parallel. Use -1 (Default) to use all but one.}

\item{click.height}{Numeric. The minimum height (in frequency bins) for a detected click to be kept.
Default is \code{10}.}
}
\value{
A tibble.
}
\description{
This function processes an sound signal to detect broadband activity by identifying 'clicks' based on
time-frame-wise (i.e., column-wise) amplitude changes in the spectrogram. It computes statistics related
to click height, variance, and centroid frequency, and can plot a spectrogram with detected clicks highlighted.
The function also classifies whether the signal contains noise or insect based on the variance and centroid
frequencies of the clicks.
}
\examples{
\dontrun{
bbai_folder(path/to/folder)
}  
}
