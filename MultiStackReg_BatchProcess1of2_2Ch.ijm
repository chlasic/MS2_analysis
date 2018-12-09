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
		for (k = 1; k < fr+1; k++) {
			for (p = 1; p < 3; p++){
				if (p==1){
					selectWindow(filenames[i]);
					aa=floor(sl/2);
			 		run("Make Substack...", "channels=1"+" slices=1-"+aa+" frames="+k);
					run("Reverse");
				}
				if (p==2){
					selectWindow(filenames[i]);
					run("Make Substack...", "channels=1"+" slices="+aa+"-"+sl+" frames="+k);
				}
			 	run("MultiStackReg", "stack_1=["+title+"-"+ p*2-1 +".tif] action_1=Align file_1=c:\\CHL\\vector.txt stack_2=None action_2=Ignore file_2=[] transformation=[Rigid Body] save");
			 	if (p==1){
			 		run("Reverse");
			 	}
			 	if (p==2){
			 		run("Slice Remover", "first=1 last=1 increment=1");
			 		}
			 	
				 // align 2nd  channels using vectors gotten from 1st channel.
			 	selectWindow(filenames[i]);
			 	if (p==1){
			 		run("Make Substack...", "channels=2"+" slices=1-"+aa+" frames="+k);
					run("Reverse");
				}
				if (p==2){
					run("Make Substack...", "channels=2"+" slices="+aa+"-"+sl+" frames="+k);
				}
			 	run("MultiStackReg", "stack_1=["+title+"-"+p*2+".tif] action_1=[Load Transformation File] file_1=c:\\CHL\\vector.txt stack_2=None action_2=Ignore file_2=[] transformation=[Rigid Body]");
				if (p==1) 	{
				 	run("Reverse");
				}
				if (p==2){
				 		run("Slice Remover", "first=1 last=1 increment=1");
				 	}
			}


			 run("Concatenate...", "  title=[" + title + "-7.tif] image1=[" + title + "-1.tif] image2=[" + title + "-3.tif] image3=[-- None --]" );
			 run("Concatenate...", "  title=[" + title + "-8.tif] image1=[" + title + "-2.tif] image2=[" + title + "-4.tif] image3=[-- None --]" );

			 	
			 // Merge into one multi-channel z-stack
			 run("Merge Channels...", "c7=["+title+"-7.tif] c6=["+title+"-8.tif] create");

			 print(k+"-th timepoint of total "+fr+" timepoints is processed ("+i+1+"th image).");
			 selectWindow("Composite");
			 	run("Arrange Channels...", "new=21");
			 if(k < 10) {
			 	selectWindow("Composite");
			 	saveAs("Tiff", dirD + "temp"+i+"/" + "fr0" + k + ".tif" );
			 	selectWindow("fr0"+k+".tif");
			 }
			 else {
			 	selectWindow("Composite");
			 	saveAs("Tiff", dirD + "temp"+i+"/" + "fr" + k + ".tif" );
			 	selectWindow("fr"+k+".tif");
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
		run("Properties...", "channels="+ch-1+" slices="+sl+" frames="+fr);

		saveAs("Tiff", dirD+title+"_reg"+".tif");
		close(); // close original file
	}
}
print("DONE!");
