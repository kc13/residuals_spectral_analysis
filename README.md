# residuals_spectral_analysis

This repository shares the Matlab code associated with the manuscript "Detecting rhythmic spiking through the power spectra of point process model residuals."

This represents the initial release of the code, sufficient for recreating the manuscript figures and tables, and the analyses and simulations that inform them. To recreate a figure or table, navigate to the corresponding subfolder under "figures" or "tables", and run the .m file beginning with the "fig" or "table" prefix.  Some subfolders include additional .m files that may be used to recreate relevant simulations or analyses. The code files are generally set up so that the necessary data files may be located automatically, assuming that that user maintains the repository's file organization, and executes the .m file from within the directory in which it is located. The only exceptions to this are the .m files for generating figures S11 and S17. In these cases, it is necessary to first run scripts to create data files that the figure-generation code requires. The data files are not supplied here due to their large size.

For the simulation code, note that information about the random number generator seeds used for the manuscript simulations is provided when available, but in some cases this may not be sufficient to exactly recreate identical simulation output.  This is due to the use of parfor loops in the simulation code.  See https://www.mathworks.com/help/parallel-computing/repeat-random-numbers-in-parfor-loops.html for more information.

The majority of the code included here is original to this repository.  An exception is the code located in the "Modified-Spline-Regression" folder.  This code was obtained from https://github.com/MehradSm/Modified-Spline-Regression.

For the remaining, original code, contributors include Karin Cox, Daisuke Kase, and Rob Turner, affiliated with the University of Pittsburgh. This research was funded in part by Aligning Science Across Parkinson’s (ASAP), through the Michael J. Fox Foundation for Parkinson’s Research (MJFF). 

Future releases will include code associated with analyses and simulations that did not directly contribute to figures or tables, and updated documentation.  Once a preprint is posted, the reference will also be provided here.


