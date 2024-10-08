% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbai_list.R
\name{nbai_list}
\alias{nbai_list}
\title{Narrow-Band Activity Index for audio.list in a List}
\usage{
nbai_list(
  audio.list,
  channel = "each",
  hpf = 0,
  freq.res = 50,
  cutoff = -60,
  activity.cutoff = 10,
  plot.binary.spec = FALSE,
  dark.plot = FALSE,
  plot.activity = FALSE,
  output.csv = "nbai_results.csv",
  n.cores = -1
)
}
\arguments{
\item{audio.list}{Character. A list containing the names of the wave audio.list to analyze.}

\item{channel}{Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.}

\item{hpf}{Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.}

\item{freq.res}{Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).}

\item{cutoff}{Numeric. Cutoff threshold defining the sounds that will be analyzed, in dBFS.}

\item{activity.cutoff}{Numeric. Cutoff percent activity. Only the frequency bands active equal or above this percentage will be considered as "active" in the active band statistics.}

\item{plot.binary.spec}{Logical. Whether to plot the binary spectrogram used for the analysis. Allowed only when Wave is mono or when one channel is selected from a stereo file.}

\item{dark.plot}{Logical. If true (default) a the binary spectrogram will have a black background.}

\item{plot.activity}{Logical. If true, a barplot depicting the percent activity in each frequency bin will be returned. Default is FALSE.}

\item{output.csv}{= Character. Name for the csv output.}

\item{n.cores}{The number of cores to use for parallel processing. Use \code{n.cores = -1} to use all but one core. Default is NULL (single-core processing).}
}
\value{
A list containing: 1) A binary spectrogram (if mono), 2) tibble with the Narrow-Band Activity Index (NBI) summary statistics, and 3) a tibble with NBI spectral, which number of rows equals the number of frequency bins in the analysis.
}
\description{
NBAI describes the relative amount of narrow-band persistent sound activity, like that of Cicadas and Orthopterans. This index can be used to evaluate insect activity and their influence on other soundscape metrics (e.g., summary acoustic indices).
}
\examples{
nbai_list(wave.list, channel = 'left', plot = TRUE)
}
