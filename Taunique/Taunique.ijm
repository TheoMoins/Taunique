run("Set Measurements...", "mean display redirect=None decimal=3");

dir = getDirectory("Choose a Directory");

testfile = dir + "test.csv";
resultfile = dir + "results.csv";

run("Camera setup", "offset=414.0 isemgain=false photons2adu=3.6 pixelsize=267.0");

run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=2.0 fitradius=3 method=[Maximum likelihood] full_image_fitting=false mfaenabled=false renderer=[Normalized Gaussian] dxforce=false magnification=5.0 dx=20.0 colorizez=false threed=false dzforce=false repaint=50");

run("Show results table", "action=changeColumnUnits column=x targetunits=nm");
run("Show results table", "action=changeColumnUnits column=sigma targetunits=nm");

Dialog.create("Choix des paramètres");
Dialog.addString("Valeur maximale de sigma [nm] :", "600");
Dialog.addString("Valeur minimale de offset [photons] :", "100");
Dialog.addString("Intensité maxiamale [photons] :", "100000");
Dialog.addString("Frame de début : ", "1");
Dialog.addString("Frame de fin : ", "2");
Dialog.addString("Distance maximale de fusion [nm] : ", "200");
Dialog.show();

s = Dialog.getString();
o = Dialog.getString();
i_max = Dialog.getString();
f1 = toString(parseFloat(Dialog.getString())-1);
f2 = toString(parseFloat(Dialog.getString())+1);
d = Dialog.getString();

run("Show results table", "action=filter formula=[ (frame > "+f1+") & (frame < "+f2+") & (offset > "+o+") & (sigma < "+s+") & (intensity < "+i_max+") ]");

run("Show results table", "action=merge zcoordweight=0.1 offframes=2 dist="+d+" framespermolecule=0");

run("Show results table", "action=changeColumnUnits column=x targetunits=px");
run("Show results table", "action=changeColumnUnits column=sigma targetunits=px");


run("Export results", "filepath="+testfile+" fileformat=[CSV (comma separated)] sigma=true intensity=false offset=false saveprotocol=false x=true y=true bkgstd=false id=false uncertainty=false frame=false detections=false");

selectWindow("Normalized Gaussian");
close();

str = File.openAsString(testfile);

lines = split(str,"\n");

roiManager("Reset");

for (i = 1; i < lengthOf(lines); i++){
	line = split(lines[i], ",");
	a = parseFloat(line[2]);
	makeOval(line[0]-a,line[1]-a, 2*a, 2*a);
	roiManager("Add");
}


roiManager("Multi Measure");

File.delete(testfile);

saveAs("Results", resultfile);

run("Scale...", "x=2 y=2 z=1.0 width=620 height=548 depth=200 interpolation=Bilinear average process create");

roiManager("Reset");

call("ij.Prefs.set","ij.Prefs.useNamesAsLabels", true);

for (i = 1; i < lengthOf(lines); i++){
	line = split(lines[i], ",");
	p0 = 2*(parseFloat(line[0])-2);
	p1 = 2*(parseFloat(line[1])-2);
	makeOval(p0,p1, 8, 8);
	ii = i-1;
	roiManager("Add");
	roiManager("Select", ii);
	roiManager("Rename", toString(i));
	roiManager("Deselect");
}
roiManager("Deselect");
run("Capture Image");
run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
close("B&C");
saveAs("PNG", dir+"screen.png");
close();

roiManager("Show All with labels");
run("Capture Image");
saveAs("PNG", dir+"screen_label.png");
close();

close();
close("ROI Manager");

close("Results");
close("Log");

Dialog.create("Done !");
Dialog.show();



