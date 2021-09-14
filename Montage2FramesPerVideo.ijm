// This script asks for a folder, finds all .TIF files it contains, 
// asks for a start frame and a stimulus frame, then pulls these
// 2 frames per video into a summary stack.
//
// D. Albrecht
// 20210608

dir = getDirectory("choose a directory of .TIF video files");
if (nImages==0){

	trackFolder=1;

	list = getFileList(dir);
	newlist = newArray(list.length);
	k=0;
	for (i=0; i < list.length; i++)  {
		if (endsWith(list[i], ".tif")) {
			newlist[k] = list[i];
			k++;
		}
	}
	
	list = Array.trim(newlist, k);

	startFrame=getNumber("Found "+list.length+" files. What initial frame # to use?", 10);
	stimFrame=getNumber("Found "+list.length+" files. What stimulus frame # to use?", 150);
	
}else{
	trackFolder=0;
	exit("Please close all images before running.");
}

if (trackFolder == 1) {

	setFont("Monospaced", 24);  // overlay settings
	setColor("white");
	x = 5; y = 30;
	
	for (i=0; i<=list.length-1; i++) { //loop through each .tif file
  	//for (i=0; i<=2; i++) {
  		
		open(dir + list[i]);  // Open the TIF
		title = getTitle();   // Get the filename

		// Copy the start slice and label it
		setSlice(startFrame);
		run("Duplicate...", "title=start.tif");
		drawString(title+" - Frame: "+startFrame, x, y, "black");

		// Copy the stimulus slice and label it
		selectWindow(title);	
		setSlice(stimFrame);
		run("Duplicate...", "title=stim.tif");
		drawString(title+" - Frame: "+stimFrame, x, y, "black");

		// Close the TIF stack
		selectWindow(title);
		close();

		// Assemble the summary stack
		if (isOpen("Summary.tif")) {
			run("Concatenate...", "title=Summary.tif image1=Summary.tif image2=start.tif image3=stim.tif");
		} else {
			run("Concatenate...", "title=Summary.tif image1=start.tif image2=stim.tif");
		}
  	}
}