// This code translates range of frames in stack to align images.

getDimensions(wi, he, ch, sl, fr);
titleT=getTitle();
rename(titleT+".tif");
titleT=getTitle();
title = substring(titleT, 0, lengthOf(titleT)-4);

waitForUser( "Pause","Locate Prev landmark and hit enter"); 
getCursorLoc(x1, y1, z1, flags1);

waitForUser( "Pause","Locate Current landmark and hit enter"); 
getCursorLoc(x2, y2, z2, flags2);

//fr1 = getNumber("Starting frame number",0);
//fr2 = getNumber("ending frame number",fr);

Stack.getPosition(chCurr, slCurr, frCurr)
fr1 = frCurr;
fr2 = frCurr+1;

xconv = x1 - x2;
yconv = y1 - y2;

selectWindow(titleT);
run("Make Substack...", "frames="+1+"-"+fr1-1);
selectWindow(titleT);
run("Make Substack...", "frames="+fr1+"-"+fr2-1);
selectWindow(titleT);
run("Make Substack...", "frames="+fr2+"-"+fr);
selectWindow(titleT);
close();

selectWindow(title+"-2.tif");
run("Translate...", "x="+ xconv +" y="+ yconv +" interpolation=None");


run("Concatenate...", "  title=[Concatenated Stacks] image1=["+title+"-1.tif] image2=["+title+"-2.tif] image3=["+title+"-3.tif] image4=[-- None --]");
selectWindow("Concatenated Stacks");
run("Properties...", "channels="+ch+" slices="+sl+" frames="+fr);

putp = random()*10;
print("DONE! (random variable: "+putp+")");



/*
StartSlice = (fr1-1)*ch*sl+1;
endSlice = fr2*ch*sl;

for (i = StartSlice; i < endSlice+1; i++) {
	setSlice(i);
	run("Translate...", "x="+ xconv +" y="+ yconv +" interpolation=None");
	print(i+"-th slice of total "+nSlices+" slices is processed.");
}
print("DONE!");
*/