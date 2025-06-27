WARNING: THIS PACKAGE IS STILL IN DEVELOPMENT!

Full documentation can be found in our Github Page: <https://panchorivasf.github.io/SoundEcology2/>

This new version of the classic ***soundecology*** package includes new acoustic indices and several additional parameters for the classic indices. Major updates include:

### Output format:

The output of each index calculation is now a long-format tibble (data frame).

### Spectrogram settings:

The first version of the package was developed when the `spectro()` function from the `seewave` package did not normalize spectrograms by default. In one of the early `seewave` updates, this behavior changed, and spectrograms were normalized by the maximum value within the recording, which is problematic when comparing metrics based on non-focal recordings (i.e., soundscape recordings). `SoundEcology2` uses non-normalized spectrograms by default and allows to modify some other parameters (e.g., frequency resolution).

### dB scale:

The original documentation stated that the amplitude threshold (or 'noise floor') for the ADI and AEI functions was in dBFS (decibels Full Scale). This was inaccurate, since the decibels were actually relative to the maximum amplitude in the recording (normalized by default). SoundEcology2 uses dBFS by default, which are computed based on the bit depth of the recording.

### Frequency resolution:

ADI and AEI functions now allow the user to select the frequency resolution via the `freq.res` parameter. Originally, this was hard-coded to 10 Hz per bin, meaning that for a recording with a 48 kHz sampling rate, the window size would be of 4800 samples. This overemphasizes the frequency resolution, blurring temporal features in the spectrogram, and artificially increasing band occupancy along the time axis. Moreover, if researchers aim to conduct standardized analyses across studies, using the same FFT window size would not yield the same resolution if the sampling frequency is not the same. In `SoundEcology2`, the window size (or window length) parameter is replaced by frequency resolution ('freq.res'), in Hz per bin.

### Diversity indices:

Acoustic diversity and evenness indices now allow users to specify how the proportions in each frequency band are computed. The original code calculated the proportion of cells above the threshold from all the cells within a frequency band. This is not a true an application of Shannon's diversity, which typically considers the proportion of individuals of a species relative to the entire population, including all the other species in the community. To address this, we have added a new option for calculating the proportion of active cells in the spectrogram matrix via the 'prop.den' (i.e., proportion denominator) argument:

-   Option '1' is the same as in the 'classic' ADI (i.e., within each frequency band).

-   Option '2' calculates the proportion over the entire range of frequencies selected by the user for analysis (i.e., min.freq to max.freq), considering all the cells above the threshold in the matrix as the denominator. This option is closer to a "real Shannon's diversity".

We also included the Frequency-dependent Acoustic Diversity Index by Xu et al. (2023), which obtains a "floating" noise profile (i.e., row-wise) for noise reduction before calculating the Acoustic Diversity Index and it doesn't use a normalized spectrogram. Alternatively it can take a noise sample to reduce noise from the analyzed files. Three functions are available to calculate this index:

-   `frequency_dependent_acoustic_diversity()`: The original implementation by Xu et al. (2023).

-   `fadi():` A harmonized implementation that yields the same output format as for the other indices in SoundEcology2.

-   `fadi_folder():` A batch option for `fadi()`.

### Channel selection:

Users can now select which channel to analyze, with an option to merge both channels into a mono file. This can potentially speed up index processing and subsequent analyses.

### New indices:

In addition to the original indices (i.e., ADI, AEI, ACI, BI, and NDSI), SoundEcology2 allows to calculate **Frequency Cover Indices (FCI)**:

*Low-Frequency Cover (LFC), Mid-Frequency Cover (MFC), High-Frequency Cover (HFC),* with user-defined frequency bands. The calculation of these indices differ from Towsey (2017) in that we don't apply noise reduction to the spectrogram, which may alter the original energy distribution of target signals. Instead, we use a binary spectrogram, where only the cells above the cutoff threshold are considered. We also included an *Ultra-Frequency Cover (UFC)* to measure the occupancy of an ultrasonic band, when available.

Three experimental acoustic indices are introduced: 1) Narrow-band Activity Index (NBAI): summarizes the percentage of active cells in each frequency bin, allowing to monitor persistent sound sources (e.g., crickets, cicadas); 2) Broad-band Activity Index (BBAI): summarizes the 'clicks' in a spectrogram, defined as clusters of energy along a column (time frame), allowing to monitor the activity of 'noisy' insects (e.g., cicadas and some katydids) and their effect on other indices, as well as geophonic noise sources such as rain and heavy wind which are generally outliers. 3) Trill Activity Index (TAI), similar to ACI, summarizes the variability in sound energy across frequency bins, being sensitive to stridulations such as those of bush-crickets and katydids. These indices can also be extracted as spectral indices (i.e., a vector of length equal to the number of frequency bins), which can be used to craft False-color Spectrograms.

-   Two function types are available for calculating acoustic indices (here 'index' replaces the name of each index):

    -   The basic function `index()` allows users to try the index with a single Wave object. It requires the file to be imported into R with the `readWave()` function from `tuneR`.
    -   The 'folder' function, `index_folder()`, takes the path to a directory and analyzes all the WAV files inside the folder. Alternatively, a list with selected files in the current working directory can be provided.

------------------------------------------------------------------------

Finally, we added several functions to facilitate common tasks, including:

`list_waves()`: Generates a list with all the WAV files in a directory.

`list_csvs()`: Generates a list with all the CSV files in a directory.

`wave_integrity()`: A function to check the integrity of WAV files in a directory. It produces a report with the last day of complete recordings and days with corrupted files, along with a corresponding plot.

`import_stereo_mix()`: Imports a stereo wave file as a mono Wave object by mixing both channels.

`import_indices()`: Imports a batch of CSV files with index data and binds them into a single tibble.

`import_kaleidoscope()`: Imports an Excel file generated with Kaleidoscope (Wildlife Acoustics). Supported formats are CSV, XLS, and XLSX.

`harmonize_index()`: Harmonizes the tibbles from `SoundEcology2` and the legacy `soundecology` package.

`var_diel_spec():` Generates a Variance Diel Spectrogram (VDS) for each day's worth of recordings in a folder. The output images are in PNG format.

`ts_diel()`: Visualizes a day's worth of acoustic indices in an interactive, smoothed time series using the Plotly package. The plots can be stored as interactive HTML files or static PNGs.

`ts_long()`: Plots a long time-series of an acoustic index, summarizing the data by day, week, and month in separate tibbles.

`ts_plus_images()`: Plots a time-series with pop-up spectrograms (e.g., generated with `var_diel_spec()`, false-color spectrograms, or other images). The path to the folder containing the spectrograms should be provided. The image files must contain the date in the format "YYYYMMDD" before the extension.

## How to install:

-   If remotes is not installed, use:

``` r
install.packages('remotes')
```

``` r
remotes::install_github("panchorivasf/SoundEcology2", dependencies = TRUE)
```
