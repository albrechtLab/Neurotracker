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
Simply copy the 'NeuroTrackerSummary.m' and 'databrowseS.m' files to a convenient MATLAB folder in your operating system.  Run the code, and point the program at the folder with the outputs of the NeuroTracker software.







