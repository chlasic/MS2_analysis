// This macro runs StackReg on all the files ending with .tif in a folder
// It seperates each channels and time points to only align one particular z-stack.
// Author: ChangHwan (Leo) Lee,     Aug 14, 2017

dirS = getDirectory("Choose source Directory");
dirD = getDirectory("Choose destination Directory");
filenames = getFileList(dirS);
for (i = 0; i < filenames.length; i++) {
	print("Working on "+i+1+"-th image of total "+filenames.length+" images");
	currFile = dirS+filenames[i];
	//dotIndex = indexOf(filenames[i], "."); 
    //title = substring(filenames[i], 0, dotIndex); 
    title = substring(filenames[i], 0, lengthOf(filenames[i])-4);
    //print(title);
     
	if(endsWith(currFile, ".tif")) { // process only files ending with .tif
		open(currFile);
		getDimensions(w,h,ch,sl,fr);
		File.makeDirectory(dirD+"temp"+i)

		// Make z-projection and align across time points
		run("Make Substack...", "channels=1");                                          // choose channel to use for alignment (after z-projection)
		selectWindow(title+"-1.tif");
		run("Z Project...", "projection=[Max Intensity] all");

//		run("Subtract...", "value=5 stack");     // added on 12/23/17 to eliminate bg glaring from oocytes

		run("MultiStackReg", "stack_1=["+"MAX_"+title+"-1.tif] action_1=Align file_1=c:\\CHL\\vector.txt stack_2=None action_2=Ignore file_2=[] transformation=[Rigid Body] save");
		selectWindow("MAX_"+title+"-1.tif");
		close();	
		selectWindow(title+"-1.tif");
		close(); 
		
		for (m = 1; m < sl+1; m++) {
			for (j = 1; j < ch+1; j++) {
				selectWindow(filenames[i]);
			 run("Make Substack...", "channels="+j+" slices="+m);
			 run("MultiStackReg", "stack_1=["+title+"-"+j+".tif] action_1=[Load Transformation File] file_1=c:\\CHL\\vector.txt stack_2=None action_2=Ignore file_2=[] transformation=[Rigid Body]");
			}
			 // Merge into one multi-channel z-stack
			 run("Merge Channels...", "c7=["+title+"-1.tif] c6=["+title+"-2.tif] create");
			 print("		"+m+"-th slice of total "+sl+" slices is processed ("+i+1+"th image).");
			 	if (m < 10) {
			 	selectWindow("Merged");
			 	run("Arrange Channels...", "new=21");
			 	saveAs("Tiff", dirD + "temp"+i+"/" + "sl0" + m + ".tif" );
			 	selectWindow("sl0"+m+".tif");
			 	}
			 	else {
			 	selectWindow("Merged");
			 	run("Arrange Channels...", "new=21");
			 	saveAs("Tiff", dirD + "temp"+i+"/" + "sl" + m + ".tif" );
			 	selectWindow("sl"+m+".tif");
			 	}
				close(); // close registered file
			}
			selectWindow(filenames[i]);
			close();

		// concatenate all image frames into one stack.
		fnames = getFileList(dirD+"temp"+i);		
		for (m = 0; m <fnames.length; m++) {
			currTfile = dirD+"temp"+i+"/"+fnames[m];
			if(endsWith(currTfile, ".tif")) { // process only files ending with .tif
				open(currTfile);
				if(m > 0) {
					run("Concatenate...", "all_open title=[Concatenated Stacks]");
				}
			}
		}

		// One big stack to hyperstack
		run("Properties...", "channels="+ch+" slices="+fr+" frames="+sl);
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		saveAs("Tiff", dirD+title+"_reg"+".tif");
		close(); // close original file

		for (m = 0; m <fnames.length; m++) {
			currTfile = dirD+"temp"+i+"/"+fnames[m];
			if(endsWith(currTfile, ".tif")) { 
				File.delete(currTfile);
			}
		}

		File.delete(dirD+"temp"+i);
		
		
	}
}
print("DONE!");
