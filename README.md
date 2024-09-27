This new version of the classic "soundecology" package includes new acoustic indices and several additional parameters for the classic indices. 
Major updates include:
-Spectrogram settings using seewave's spectro() function are now fully customizable, which has several implications. 
The first version of the package was developed when the spectro() function from the seewave package did not normalize the spectrograms by default.
In one of the early seewave updates, this behavior changed and then the spectrograms were normalized by the maximum value within the recording, which 
is problematic if you want to compare metrics based on non-focalized recordings (i.e., soundscape recordings). Additionally, the original documentation 
states that the decibel threshold (or 'noise floor') for functions such as ADI and AEI was measured in dBFS (decibels Full Scale), which was inacurate, since 
the decibels were relative to the maximum amplitude in the recording (normalized). In SoundEcology2 the user can select dBFS when using non-normalized spectrograms 
(which is highly recommended if the goal is to obtain valid comparisons across recordings). 
Additionally, the ADI and AEI functions now allow the user to select the frequency resolution ('freq.res' parameter). Originally, this was hard-coded to 10 Hz per bin, 
meaning that for a recording with a sampling rate of 48 kHz, the window size will be 4800 samples, which overemphasizes the frequency resolution, blurring time features in the spectrogram. 
These diversity indices now allow the user to choose between diffrent ways of calculating the proportions in each frequency band. The original code takes the proportion of cells above the threshold
within the context of each frequency band. Some argue that this is not trully an application of Shannon's diversity, since the classical equation takes the proportion
of individuals in the whole population, including all the species in the community. Therefore, we have added two new options for calculating proportion of active cells 
in the spectrogram matrix: Option '2' calculates the proportion over the whole range of frequencies selected by the user for the analysis, which would be closer to a 
"real Shannon's Diversity", while Option '3' uses all the frequencies in the spectrogram (up to the Nyquist frequency), which could be useful when the goal is to detect
unusual activity like heavy rain. Another issue we wanted to address is to find a way to produce spectrograms with a standard resolution, regardless of the sampling rate. 
If researchers want to conduct standardized analyses across studies, using the same window size (or window length) would not return the same resolution if the sampling 
frequency was different. Therefore, we decided to provide always the frequency resolution as a parameter instead of window size. 

-The output of each index calculation is now a  long-format tibble (data frame).
-The user can select wich channel to analyze, including an option to parse a 'mix' of the two channels into a mono file, potentially increasing processing and posterior analysis efficiency. 
-Each index comes with 3 function types to make things easier:
1) The basic function "index()" allows you to try the index with a single file. 
2) The 'list' function type "index_list()" receives a list with the names of WAV files in the working directory, so you can experiment with small batches of files. 
3) The 'folder' function type "index_folder()" takes the path to a directory and analyzes all the WAV files inside the folder. 

Finally, we added some helper functions to help speed up the analyses, including:
list_waves : a function to list all the wave files in a directory.
wave_integrity: a function to check the integrity of wave files in a directory. It produces a report with the last day with complete recordings and the days with corrupted files, along with a plot with these data.
diel_index_plot: a function to visualize a day's worth of acoustic indices in an interactive smoothed (LOESS) time series using the Plotly package. The plots can be stored as an interactive HTML file or static PNG.
