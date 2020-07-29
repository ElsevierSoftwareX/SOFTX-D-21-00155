r"""
sas_temper

test_datareading.py:  code for testing reading of a data file  

Read in a 1D sans data set and confirm that the code works by outputting plots.

This is intended for the "tests" folder.  

Oak Ridge National Laboratory, 2020

"""

from matplotlib import pyplot as plt
import sas_temper.sas_data as sas_data
import numpy as np
import sys

def main():
	if len(sys.argv) < 2:
		raise Exception("No configuration file specified")
	
	plot_model(sys.argv[1], showPlot=False, saveImage=True)
	
def plot_model(filename, showPlot = False, saveImage = True):
	#this clearly does not read a data file at this point in time
	m1 = sas_data.SASData(filename,qmin=0.004,qmax=0.500)
	m2 = sas_data.SASData(filename,0.010,0.250)
	m3 = sas_data.SASData(filename,0.040,0.150)
	
	#do some scaling to offset the curves
	m2.y *= 10.0
	m3.y *= 100.0
	
	#make a figure of the data and save it
	fig = plt.figure(figsize = [4,3],dpi = 100)
		
	grph = fig.add_subplot(1,1,1)
	grph.loglog(m1.x,m1.y,'bo', m2.x,m2.y,'ro', m3.x,m3.y,'go', nonposy = 'mask')
	grph.set_title( "test_datareading.py" )
	
	
	# show the image, if desired
	if showPlot :
		plt.show()
	
	# save the image, if desired
	if saveImage :
		oname = str(filename) + ".png"
		fig.savefig(oname,format="png")

		
if __name__=="__main__":
	main()
