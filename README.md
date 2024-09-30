This new version of the classic **"soundecology"** package includes new acoustic indices and several additional parameters for the classic indices.
Major updates include:

### Spectrogram settings: 

Using the seewave package's spectro() function, the spectrogram settings are now fully customizable. The first version of the package was developed when the spectro() function from the seewave package did not normalize spectrograms by default. In one of the early seewave updates, this behavior changed, and spectrograms began to be normalized by the maximum value within the recording, which is problematic when comparing metrics based on non-focal recordings (i.e., soundscape recordings). SoundEcology2 uses non-normalized spectrograms by default. 

### dB scale: 

The original documentation stated that the decibel threshold (or 'noise floor') for functions such as ADI and AEI was measured in dBFS (decibels Full Scale), which was inaccurate, as the decibels were actually relative to the maximum amplitude in the recording (normalized). Additionally, if one wants to use non-normalized values in a dB scale, the range will include both negative and positive values, whwich is counterintuitive when one wants to select a dB cutoff threshold. SoundEcology2, uses dBFS by defaultm, which is calculated based on the bit depth of the recordings.


### Frequency resolution: 

ADI and AEI functions now allow the user to select the frequency resolution via the freq.res parameter. Originally, this was hard-coded to 10 Hz per bin, meaning that for a recording with a 48 kHz sampling rate, the window size would be 4800 samples. This overemphasized the frequency resolution, blurring temporal features in the spectrogram, potentially increasing occupancy along the time axis as a side effect. Moreover, if researchers aim to conduct standardized analyses across studies, using the same window size (or window length) would not yield the same resolution if the sampling frequency differed. To produce spectrograms with a standard resolution regardless of the sampling rate, we can decouple window size from frequency resolution. In SoundEcology2, the window size (or window length) parameter is replaced by frequency resolution ('freq.res'), in Hz per bin. 

### Diversity indices: 
These indices now allow users to specify how the proportions in each frequency band are computed. The original code calculated the proportion of cells above the threshold within each frequency band. Some argue that this is not truly an application of Shannon's diversity, which typically considers the proportion of individuals of a species relative to the entire population, including all the other species in the community. To address this, we have added two new options for calculating the proportion of active cells in the spectrogram matrix via the 'prop.den' (proportion denominator) argument:


- Option '1' is the same as in the 'classic' ADI (i.e., within each frequency band). 

- Option '2' calculates the proportion over the entire range of frequencies selected by the user for analysis, which is closer to "real Shannon's diversity".


- Option '3' uses all the frequencies in the spectrogram (up to the Nyquist frequency), which can be useful for detecting unusual activities, such as heavy rain.
Output format: The output of each index calculation is now a long-format tibble (data frame).

### Channel selection: 

Users can now select which channel to analyze, with an option to merge both channels into a mono file. This can potentially speed up index processing and subsequent analyses.

### New indices: 

In addition to the original indices (i.e., ADI, AEI, ACI, BI, and NDSI), SoundEcology2 allows to calculate **Frequency Cover Indices (FCI)**: 

*Low-Frequency Cover (LFC), Mid-Frequency Cover (MFC), High-Frequency Cover (HFC), and Ultra-Frequency Cover (HFC)* with user-defined frequency bands. The calculation of these indices differ from Towsey (2017) in that we don't apply noise reduction to the spectrogram, which may alter the original energy distribution of target signals. Instead, we use a binary spectrogram, where only the cells above the cutoff threshold are considered. Three unpublished acoustic indices are introduced: 1) Narrow-band Activity Index (NBAI): summarizes the percentage of active cells in each frequency bin, allowing to monitor persistent sound sources (e.g., crickets, cicadas); 2) Broad-band Activity Index (BBAI): summarizes the active cells along each time frame, allowing to monitor the activity of 'noisy' insects (e.g., cicadas and some katydids) and their effect on other indices, as well as geophonic noise sources such as rain and heavy wind which are generally outliers. 3) Trill Activity Index (TAI), similar to ACI, summarizes the variability in sound energy across frequency bins, being more sensitive to stridulations such as those of bush-crickets and katydids. These indices can also be extracted as spectral indices (i.e., a vector of length = number of frequency bins), which can be used to craft False-color Spectrograms. 

- Three function types for each index:

    - The basic function index() allows users to try the index with a single file.

    - The 'list' function, index_list(), accepts a list of WAV file names in the working directory, enabling users to experiment with small batches of files.

    - The 'folder' function, index_folder(), takes the path to a directory and analyzes all the WAV files inside the folder.

--- 
Finally, we have added some helper functions to streamline the analysis process, including:

list_waves: A function to list all the WAV files in a directory.

wave_integrity: A function to check the integrity of WAV files in a directory. It produces a report with the last day of complete recordings and days with corrupted files, along with a corresponding plot.

diel_index_plot: A function to visualize a day's worth of acoustic indices in an interactive, smoothed (LOESS) time series using the Plotly package. The plots can be stored as interactive HTML files or static PNGs.

## How to install:

- If remotes is not installed, use: install.packages('remotes').

- library(remotes)

- remotes::install_github("panchorivasf/SoundEcology2", dependencies = TRUE)
