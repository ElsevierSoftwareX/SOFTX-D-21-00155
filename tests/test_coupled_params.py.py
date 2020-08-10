r"""
sas_temper

test_coupled_params.py:  code for testing parameters that are coupled  

This is intended for the "tests" folder.  

Oak Ridge National Laboratory, 2020

"""

import sys
import sas_temper.sa_config as sa_config
import sas_temper.parse_conf as parse_conf
import sas_temper.modelconfig as modelconfig
import sas_temper.param as param
import sas_temper.polydispersity as polydispersity
import sas_temper.sas_temper_engine as engine

def main():
	if len(sys.argv) < 2:
		raise Exception("No configuration file specified")
	
	# this now works like the main function in __main__.py
	sasTemperConf, modelConf = parse_conf.parse_config(sys.argv[1])
	
	print("the input data file is " + str(sasTemperConf.datafile))
	print("    to be fit over q = " + str(sasTemperConf.qmin) + " to " + str(sasTemperConf.qmax))
	print("the output filename base is " + str(sasTemperConf.output))
	print("")
	print("The parameters for the simulated annealing:")
	print("     temperatures: " + str(sasTemperConf.temperatures))
	print("     temperature_rate: " + str(sasTemperConf.temp_rate))
	print("     parameter_rate: " + str(sasTemperConf.param_rate))
	print("     iterations: " + str(sasTemperConf.iterations))
	print("     models to generate: " + str(sasTemperConf.models))
	
	print("")
	print("The model that will be used for the data:")
	print("     name:  " + str(modelConf.name))
	print("     category:  " + str(modelConf.category))
	
	for i,p in enumerate(modelConf.params):
		print("\t\tparameters name:  " + str(p.name) + " kind: " + str(p.kind))
		if p.kind in ["fixed"]:
			print("\t\t\t value:  " + str(p.min))
		else :
			(print "\t\t\t range:  " + str(p.min) + " to " + str(p.max))
		if p.coupled is not None :
			print("\t\t\t coupled to " + str(p.coupled))
	
	if modelConf.sq is not None:
		print("\tstructure factor:  " + modelConf.sq.type)
		for i,sqp in enumerate(modelConf.sq.params):
			print("\t\tparameters name:  " + str(sqp.name) + " kind: " + str(sqp.kind))
			if sqp.kind in ["fixed"]:
				print("\t\t\t value:  " + str(sqp.min))
			else :
				print("\t\t\t range:  " + str(sqp.min) + " to " + str(sqp.max))
			if sqp.coupled is not None :
				print("\t\t\t coupled to " + str(sqp.coupled))
	
    # create a couple instances of the model configuration class to hold things
    f = modelconfig.ModelConfig(modelConf.name,modelConf.category,modelConf.params,modelConf.sq)
    conv_conf = modelconfig.ModelConfig(modelConf.name,modelConf.category,modelConf.params,modelConf.sq)
    
    # then call the function that will create a single random model
    f = engine.define_model(1,modconf,10.0,1.0,cur)
    
    # convert it to a configuration that is used for calculating an actual model profile
    conv_conf = modelconfig.convert_conf(f)
    
    # print out the results
    for i,p in enumerate(modelConf.params):
		print("\t\tparameters name:  " + str(p.name) + " value: " + str(p.val))
	
	if modelConf.sq is not None:
		print("\tstructure factor:  " + modelConf.sq.type)
		for i,sqp in enumerate(modelConf.sq.params):
			print("\t\tparameters name:  " + str(sqp.name) + " value: " + str(sqp.val))
			
if __name__=="__main__":
	main()
