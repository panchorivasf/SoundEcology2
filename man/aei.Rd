% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aei.R
\name{aei}
\alias{aei}
\title{Acoustic Evenness Index}
\usage{
aei(
  wave,
  freq.res = 50,
  w.fun = "hanning",
  min.freq = 0,
  max.freq = 10000,
  n.bands = 10,
  cutoff = -60,
  norm.spec = FALSE,
  noise.red = 0,
  rm.offset = TRUE,
  props = TRUE,
  prop.den = 1,
  db.fs = TRUE
)
}
\arguments{
\item{wave}{an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.}

\item{w.fun}{window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".}

\item{min.freq}{minimum frequency to compute the spectrogram.}

\item{max.freq}{maximum frequency to compute the spectrogram.}

\item{n.bands}{number of bands to split the spectrogram.}

\item{cutoff}{dB threshold to calculate energy proportions.}

\item{norm.spec}{logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended because normalized spectrograms with different SNR are not comparable).}

\item{noise.red}{numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.}

\item{rm.offset}{logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.}

\item{props}{logical; if set to TRUE, the function stores the energy proportion values for each frequency band and channel. Default = TRUE.}

\item{prop.den}{numeric; indicates how the energy proportion is calculated.}

\item{db.fs}{logical; if TRUE, the amplitude scale is expressed as decibels Full Scale (dBFS). Only used when norm = FALSE.}

\item{frew.res}{the frequency resolution  (Hz per bin) to use. From this value the window length for the FFT will be calculated (sampling rate / frequency resolution).}
}
\value{
A tibble (data frame) with the AEI values for each channel (if stereo), metadata, and the parameters used for the calculation.
}
\description{
Acoustic Evenness Index from Villanueva-Rivera \emph{et al.} 2011.
The AEI is calculated by dividing the spectrogram into frequency bands (default 10),
taking the proportion of energy in each band above an energy threshold,
and then calculating the Gini Coefficient from those proportions.
The new version allows the user to choose between different ways to compute
the proportions before calculating the Gini, among other new parameter options (see Details)
}
\examples{
aei(tropicalsound)
}
