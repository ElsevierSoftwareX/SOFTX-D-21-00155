r"""
sas_temper

test_smearing.py:  perform a test calculation that reads in real
					data with smearing and then calculates the 
					smeared and unsmeared profiles and make a graph
					of both of with the input data.

This is intended for the "tests" folder.  

Oak Ridge National Laboratory
William Heller, 2017

"""

import sys
import numpy as np
from sasmodels.sasview_model import _make_standard_model

import sas_temper.sas_data as data
import sas_temper.sas_calc as calc
import sas_temper.parse_conf as parse_conf
import sas_temper.modelconfig as modelconfig
import sas_temper.sas_temper_config as config

from matplotlib import pyplot as plt

def main():
	if len(sys.argv) < 2:
		raise Exception("No configuration file was specified")
	
	# parse the configuration file for the parameters
	#modelConf = modelconfig.ModelConfig()
	#sasTemperConf = sa_config.SAConfiguration()
	sasTemperConf, modelConf = parse_conf.parse_config(sys.argv[1])
	
	# get the input data specified in the configuration file
	# this also creates the model intensity profiles
	input_data = data.SASData(sasTemperConf.datafile,np.float64(sasTemperConf.qmin),np.float64(sasTemperConf.qmax))
	smeared_calc = data.Model(input_data, unsmeared=False)
	unsmeared_calc = data.Model(input_data, unsmeared=True)
	
	# the parameters for the model have to be set here
	# for this test, parameter setting is going to be ugly
	for i,p in enumerate(modelConf.params):
		if p.name is "scale":
			modelConf.params[i].val = 1.00
		if p.name is "background":
			modelConf.params[i].val = 0.001
		if p.name is "sld":
			modelConf.params[i].val = 8.00
		if p.name is "sld_solvent":
			modelConf.params[i].val = -0.56
		if p.name is "radius":
			modelConf.params[i].val = 100.0
			
	# now, we do the calculation
	print "model configuration name " + str(modelConf.name)
	print "input data file "+str(sasTemperConf.datafile) + " q-range "+str(sasTemperConf.qmin)+"\t"+str(sasTemperConf.qmax)
	print "\t\tinput data dq "+str(input_data.dx[10])+" at a q of "+str(input_data.x[10])
	print "\t\tlength of the input data smearing array "+str(len(input_data.dx))
	
	unsmeared_calc = calc.calc_profile_usm(input_data, modelConf)
	smeared_calc = calc.calc_profile(input_data, modelConf, unsmeared_calc)
		
	# now, make a figure of the three profiles and save it
	fig = plt.figure(figsize = [4,3],dpi = 200)
		
	grph = fig.add_subplot(1,1,1)
	grph.loglog(input_data.x,input_data.y,'b-', smeared_calc.x,smeared_calc.y,'r-', unsmeared_calc.x,unsmeared_calc.y,'g-', nonposy = 'mask')
	grph.set_title( "test_smearing.py" )
	oname = "test_smearing.png"
	fig.savefig(oname,format="png")

	
if __name__=="__main__":
	main()
