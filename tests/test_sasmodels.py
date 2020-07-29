r"""
sas_temper

test_sasmodels.py:  perform a test calculation that calls SASModels
				and create a plot of the results for easy review

This is intended for the "tests" folder.  It is only intended
to make sure that SASModels is around, that it is being called correctly,
and that the code is working reasonably correctly.  

Oak Ridge National Laboratory, 2020

"""

import sys
import numpy as np
from sasmodels.sasview_model import _make_standard_model

import sas_temper.sas_data as data

from matplotlib import pyplot as plt

def main():
	if len(sys.argv) < 2:
		raise Exception("No input sans data file specified")
		
	input_data1 = data.SASData(sys.argv[1],qmin=0.004,qmax=0.50)
	input_data2 = data.SASata(sys.argv[1],qmin=0.01,qmax=0.30)
	input_data3 = data.SASData(sys.argv[1],qmin=0.04,qmax=0.15)
	output_data1 = data.Model(input_data1, unsmeared=False)
	output_data2 = data.Model(input_data2, unsmeared=True)
	output_data3 = data.Model(input_data3, unsmeared=False)

	Model = _make_standard_model('sphere')
	model_sphere = Model()
	sld = 6.38
	sld_solvent = -0.56
	radius = 100.0
	scale = 1.0
	background = 0.001
	model_sphere.setParam('sld',sld)
	model_sphere.setParam('sld_solvent',sld_solvent)
	model_sphere.setParam('radius',radius)
	model_sphere.setParam('scale',scale)
	model_sphere.setParam('background',background)
	output_data1.y = model_sphere.evalDistribution(output_data1.x)
	
	#model_sphere.setParam('sld',sld)
	#model_sphere.setParam('sld_solvent',sld_solvent)
	#model_sphere.setParam('radius',radius)
	#model_sphere.setParam('scale',scale)
	#model_sphere.setParam('background',background)
	output_data2.y = model_sphere.evalDistribution(output_data2.x)
	output_data2.y *= 10.0
	
	#model_sphere.setParam('sld',sld)
	#model_sphere.setParam('sld_solvent',sld_solvent)
	#model_sphere.setParam('radius',radius)
	#model_sphere.setParam('scale',scale)
	#model_sphere.setParam('background',background)
	output_data3.y = model_sphere.evalDistribution(output_data3.x)
	output_data3.y *= 100.0
	
	#make a figure of the data and save it
	fig = plt.figure(figsize = [4,3],dpi = 100)
		
	grph = fig.add_subplot(1,1,1)
	grph.loglog(output_data1.x,output_data1.y,'bo', output_data2.x,output_data2.y,'ro', output_data3.x,output_data3.y,'go', nonposy = 'mask')
	grph.set_title( "test_sasmodels.py" )
	oname = "test_sasmodels.png"
	fig.savefig(oname,format="png")

	
if __name__=="__main__":
	main()
