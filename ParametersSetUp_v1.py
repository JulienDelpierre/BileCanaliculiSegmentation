# @ String (visibility = MESSAGE, value="<html>##########################################################################################<br/>__________________________<b> File to test </b>________________________________<br/>##########################################################################################</html>") docmsg1
# @ File    (label = "Input File for testing", style = "file") srcFile
# @ Boolean (label = "Open with bioformat" , value=True) bifoboo
# @ Boolean (label = "Convert stack to hyperstack" , value=True) convboo
# @ Integer (label = "Nb of Channels (if convertion is needed)", value="4") nch
# @ Integer (label = "Nb of Z Slices (if convertion is needed)", value="16") nsl
# @ String (visibility = MESSAGE, value="<html>##########################################################################################<br/>________________________________<b> Segmentation parameters </b>___________________________________<br/>##########################################################################################</html>") docmsg2
# @ String (visibility = MESSAGE, value="<html><b>Remove saturating pixels</b></html>") docmsga
# @ String  (label = "List of channels to clean (format: 2,3 / to skip: 0)", value="2,3") dirtlist
# @ Integer (label = "Maximum intensity (for cleaning)", style="slider", min=0, max=65500, stepSize=1000, value=50000) dirtthre
# @ String (visibility = MESSAGE, value="<html><b>Apical Surfaces</b></html>") docmsgz
# @ Integer (label = "Apical Surface Marker Channel (0 if none)", value="2") chcd
# @ Integer (label = "Apical Surface signal smoothing", style="slider", min=0, max=10, stepSize=1, value=2) cdsmo
# @ String (visibility = MESSAGE, value="<html><b>Tight Junctions</b></html>") docmsge
# @ Integer (label = "Junctions Marker Channel (0 if none)", value="3") chju
# @ Integer (label = "Junctions Marker signal smoothing", style="slider", min=0, max=10, stepSize=1, value=2) jusmo
# @ String (visibility = MESSAGE, value="<html><b>Actin</b></html>") docmsgr
# @ Integer (label = "Phalloidin Channel (0 if none)", value="1") chph
# @ Integer (label = "Phalloidin signal smoothing", style="slider", min=0, max=10, stepSize=1, value=5) phsmo

import os
from ij import IJ, ImagePlus, ImageStack
from ij.plugin import ImageCalculator
from ij import WindowManager
from ij import process as proc

ic = ImageCalculator()
w = WindowManager

####Interpret user input
dirtboo = True
if dirtlist == '0':
	dirtboo = False
dirtlist = list(dirtlist.split(','))
dirtlist = [int(a) for a in dirtlist]

cdboo=True
if chcd == 0:
	cdboo=False

juboo=True
if chju == 0:
	juboo=False

phboo=True
if chph == 0:
	phboo=False
####

def run():
	if bifoboo:
		IJ.run("Bio-Formats Importer", "open="+str(srcFile)+" autoscale color_mode=Default view=Hyperstack stack_order=XYCZT")
		imp = IJ.getImage()
		imp.setTitle("OriginalFile")
		process("OriginalFile", nch, nsl)
	if not bifoboo:
		IJ.open(str(srcFile))
		imp = IJ.getImage()
		imp.setTitle("OriginalFile")
		process("OriginalFile", nch, nsl)

###########################Helper Functions
def sourcimage(nch, nsl):
	imp = IJ.getImage()
	if convboo:
		nsli = imp.getNSlices()
		nfra = nsl
		if nsl != nsli and nsl>1 :
			IJ.run("Stack to Hyperstack...","order=xyczt(default) channels="+str(nch)+" slices="+str(nsl)+" frames=1 display=Color")
		imp = IJ.getImage()
	if not convboo:
		nfra = imp.getNSlices()
		nch = imp.getNChannels()
	calib = imp.getCalibration()
	calib.setUnit("um")
	nfra = imp.getDimensions()[3]
	imp = IJ.getImage()
	imp.setTitle("source")

#Get Maximal pixel value and associated slice from a stack
def getmxint(imp, nsl):
	slt = 0
	intt = 0
	for i in range(1,nsl+1):
		IJ.setSlice(i)
		ip = imp.getProcessor().convertToFloat()  
		pixels = ip.getPixels() 
		mx = max(pixels)
		if mx > intt:
			intt = mx
			slt = i
	return intt, slt

#Remove saturating pixels and surrounding area
def RemoveDirt(cibl, cha, sortie, nch, nsl):
	IJ.selectWindow(cibl)
	if dirtboo:
		nampara = ''
		for q in range(nch):
			l = q+1
			nampara = nampara+'c'+str(l)+'=Masked'+str(l)+' '
			if l in cha:
				IJ.selectWindow(cibl)
				IJ.run("Duplicate...", "duplicate channels="+str(l))
				imp = IJ.getImage()
				imp.setTitle("Source"+str(l))
				imp = IJ.getImage()
				IJ.run("Duplicate...", "duplicate channels="+str(l))
				imp = IJ.getImage()
				imp.setTitle("Msk")
				imp = IJ.getImage()
				intt, slt = getmxint(imp, nsl)
				if intt < dirtthre:
					imp.setTitle("Masked"+str(l))
					IJ.selectWindow("Source"+str(l))
					clos()
					continue
				IJ.setSlice(slt)
				IJ.setThreshold(dirtthre, 65535)
				IJ.run("Convert to Mask", "method=Default background=Dark")
				for a in range(5):
					IJ.run("Dilate (3D)", "iso=255")
				IJ.setThreshold(0, 1)
				IJ.run("Convert to Mask", "method=Default background=Dark")
				imp = IJ.getImage()
				IJ.run(imp, "Divide...", "value=255.000 stack")
				imstk2 = w.getImage("Source"+str(l))
				imstk1 = w.getImage("Msk")
				imp3 = ic.run("multiply stack", imstk2, imstk1)
				IJ.selectWindow("Source"+str(l))
				imp = IJ.getImage()
				imp.setTitle("Masked"+str(l))
				IJ.selectWindow("Msk")
				clos()
			else:
				IJ.selectWindow(cibl)
				IJ.run("Duplicate...", "duplicate channels="+str(l))
				imp = IJ.getImage()
				imp.setTitle("Masked"+str(l))

		IJ.run("Merge Channels...", nampara+"create")
		imp = IJ.getImage()
		imp.setTitle(sortie)
		
	if not dirtboo:
		IJ.selectWindow(cibl)
		imp = IJ.getImage()
		imp.setTitle(sortie)

#Close modified images without saving
def clos():
	imp = IJ.getImage()
	imp.changes = False
	imp.close()

#Basic segmentation
def SimpleSeg(cibl, chan, rad, sortie, method, slitest):
	IJ.selectWindow(cibl)
	IJ.run("Duplicate...", "duplicate channels="+chan+slitest)
	IJ.run ("Median...", "radius="+rad+" stack")
	imp = IJ.getImage()
	imp.setTitle(sortie)
	IJ.run("Auto Threshold", "method="+method+" ignore_black ignore_white white stack use_stack_histogram")

#Remove background based on gaussian blur
def pretreat(cibl, chan, sortie):
#background
	bck = "BckMask"+str(chan)
	fore = "ForMask"+str(chan)
	IJ.selectWindow(cibl)
	IJ.run("Duplicate...", "duplicate channels="+chan)
	imp = IJ.getImage()
	imp.setTitle(bck)
	IJ.run("Gaussian Blur...", "sigma=4 scaled stack"),
#Foreground
	IJ.selectWindow(cibl)
	IJ.run("Duplicate...", "duplicate channels="+chan)
	imp = IJ.getImage()
	imp.setTitle(fore)
	IJ.run ("Median...", "radius=3 stack")
#Pre-treatment
	imstk2 = w.getImage(bck)
	imstk1 = w.getImage(fore)
	imp3 = ic.run("substract stack", imstk1, imstk2)
	imp = IJ.getImage()
	imp.setTitle(sortie)
	IJ.selectWindow(bck)
	clos()

##### Main
def process(fil, nch, nsl):

	sourcimage(nch, nsl)
	IJ.selectWindow('source')
	imp = IJ.getImage()
	nch = imp.getNChannels()
	nsl = imp.getNSlices()
	IJ.run("Duplicate...", "duplicate")

	#Remove saturating pixels
	RemoveDirt('source', dirtlist, "clean", nch, nsl)

	#Segmentation on CD13 and Phalloidin
	IJ.selectWindow('clean')
	imp = IJ.getImage()
	nsl = imp.getNSlices()
	nslt = int(nsl/2)
	nsltest = " slices="+str(nslt)
	nsltest2 = " range="+str(nslt)+'-'+str(nslt)

	texttoprint='\n\nTesting the file: '+str(srcFile)+'\nList of operations:\n  Dirt Removal:  '''+str(dirtboo)+'\n  Apical surface segmentation:  '+str(cdboo)+'\n  Phalloidin segmentation:  '+str(phboo)+'\n  Junctions segmentation:  '+str(juboo)+'''
	\nParameters:\n  List of channels to clean:  '''+str(dirtlist)+'\n  Upper threshold for cleaning:  '+str(dirtthre)+'\n\n  Slice used for testing:  '+str(nslt)+'\nApical surface marker channel:  '+str(chcd)+'\n  Apical Surface signal smoothing:  '+str(cdsmo)+'''\nJunctions marker channel:  '''+str(chju)+'\n  Junctions marker signal smoothing:  '+str(jusmo)+'''\nPhalloidin channel:  '''+str(chph)+'\n  Phalloidin signal smoothing:  '+str(phsmo)
	
	if cdboo:
		SimpleSeg("clean", str(chcd), str(cdsmo), "Apicalsmooth", "[Try all]", nsltest)
		imp = IJ.getImage()
		imp.setTitle("ApicalSegmentationTest")
	if phboo:
		SimpleSeg("clean", str(chph), str(phsmo), "Phallosmooth", "[Try all]", nsltest)
		imp = IJ.getImage()
		imp.setTitle("PhalloidinSegmentationTest")
	if juboo:
		pretreat("clean", str(chju), "PreTreatedJunctions")
		SimpleSeg("PreTreatedJunctions", '1', str(jusmo), "Junctionsmooth", "[Try all]", nsltest2)
		imp = IJ.getImage()
		imp.setTitle("JunctionSegmentationTest")
	print(texttoprint)

run()