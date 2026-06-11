# NeuroTracker
NeuroTracker is a software plugin for Fiji/ImageJ that allows the analysis of epifluorescent imaging data for multi-slice TIFF image stacks.

# Usage
The uploaded PDF file gives a step-by-step guide to installing and running the plugin using Fiji.

General Instructions:
* Download files from this Github repo
* Place the jar file included in this archive into the plugins folder within the installation location of Fiji/ImageJ.
* Place "Montage2FramesPerVideo.ijm" and "StartupMacros.fiji.ijm" into the 'macros' subfolder within the main Fiji plugins folder.
* Start/Restart Fiji/ImageJ
* Run Neurotracker from Plugins > Tracking > Neurotracker
* Follow the instructions in the included PDF

# NeuroTrackerSummary
NeuroTrackerSummary is a MATLAB file that reads and formats the data created as outputs by the included NeuroTracker software plugin.
Copy the 'NeuroTrackerSummary.m' and 'databrowseS.m' files to a convenient MATLAB folder in your operating system.  Run the code, and point the program at the folder with the outputs of the NeuroTracker software.

For more information, see https://albrechtlab.github.io/neurotracker/

# Example video and analysis
Two neural recording trials from the publication [White et al., JoVE, 2023](https://www.jove.com/v/65042/automated-multimodal-stimulation-simultaneous-neuronal-recording-from) are available [here](https://users.wpi.edu/~dalbrecht/neurotracker/training) (warning: 1 Gb total!)

These represent the first two trials contributing to Fig. 8c, top panel, of the publication. Here, a set of >20 animals are stimulated optically and chemically withing the same 45-s trial.  
<img width="2838" height="1161" alt="65042fig08large" src="https://github.com/user-attachments/assets/90ade978-ba4d-4006-96bf-c386fe29bf7c" />

After running NeuroTracker, example track files are cintained within the **exampleTrackingFiles.zip** file in this repository.



