r"""
sas_temper

test_parse.py:  code for testing YAML parsing of config file  

Read in a YAML 1.2 compliant configuration file for the sas_temper
fitting engine that is built on simulated annealing. 

This is intended for the "tests" folder.  

Oak Ridge National Laboratory, 2020

"""

import sys
import sas_temper.sa_config as config
import sas_temper.param as param
import sas_temper.polydispersity as polydispersity
#import structurefactor

def main():
	if len(sys.argv) < 2:
		raise Exception("No configuration file specified")
	
	conf_params = config.load_config(config_file=sys.argv[1])
	
	print "the file is " + str(conf_params["file"]["name"])
	print "    to be fit over q = " + str(conf_params['file']['qmin']) + " to " + str(conf_params['file']['qmax'])
	print "the output filename base is " + str(conf_params["output_files"])
	
	sa_parameters = {}
	sa_parameters["name"] = str(conf_params["file"]["name"])
	sa_parameters["qmin"] = float(conf_params['file']['qmin'])
	sa_parameters["qmax"] = float(conf_params['file']['qmax'])
	sa_parameters["output"] = str(conf_params["output_files"])
	
	print ""
	print "The parameters for the simulated annealing:"
	print "     temperatures: " + str(conf_params["sa_parameters"]["temperatures"])
	print "     temperature_rate: " + str(conf_params["sa_parameters"]["temperature_rate"])
	print "     parameter_rate: " + str(conf_params["sa_parameters"]["parameter_rate"])
	print "     iterations: " + str(conf_params["sa_parameters"]["iterations"])
	print "     models to generate: " + str(conf_params["sa_parameters"]["models_to_generate"])
	
	sa_parameters["temperatures"] = str(conf_params["sa_parameters"]["temperatures"])
	sa_parameters["temperature_rate"] = str(conf_params["sa_parameters"]["temperature_rate"])
	sa_parameters["parameter_rate"] = str(conf_params["sa_parameters"]["parameter_rate"])
	sa_parameters["iterations"] = str(conf_params["sa_parameters"]["iterations"])
	sa_parameters["models"] = str(conf_params["sa_parameters"]["models_to_generate"])
	
	sasTemperConf = config.SAConfiguration(sa_parameters)
	
	print ""
	print "The model that will be used for the data:"
	print "     name:  " + str(conf_params["model"]["name"])
	print "     category:  " + str(conf_params["model"]["category"])
	
	model_params = []
	for i, params in enumerate(conf_params["model"]):
		if params not in ["name","category","Structure_Factor"]:
			txt1 = str(params)
			print ""
			print "     parameter name:  " + txt1
			
			coupled_name=None
			pd_name=None
			var_kind=None
			for j, mod in enumerate(conf_params["model"][txt1]):
				txt2 = str(mod)
				if mod in ["fixed", "linear", "log", "integer"]:
					var_kind = txt2
					if txt2 == "fixed":
						print "          type: " + txt2 + "; value: " + str(conf_params["model"][txt1][txt2][0])
						min_val = float(conf_params["model"][txt1][txt2][0])
						max_val = float(conf_params["model"][txt1][txt2][0])
					elif txt2 in ["linear", "log", "integer"]:
						print "          type: " + txt2 + "; range: " + str(conf_params["model"][txt1][txt2][0]) + " " + str(conf_params["model"][txt1][txt2][1])
						min_val = float(conf_params["model"][txt1][txt2][0])
						max_val = float(conf_params["model"][txt1][txt2][1])
				elif mod == "coupled":
					coupled_name = str(conf_params["model"][txt1][txt2][0])
					print "          coupled to: " + str(conf_params["model"][txt1][txt2][0])
				elif mod == "polydispersity":
					for k, poly in enumerate(conf_params["model"][txt1][txt2]):
						txt3 = str(poly)
						pd_name = txt3
						pd_minval = float(conf_params["model"][txt1][txt2][txt3][0])
						pd_maxval = float (conf_params["model"][txt1][txt2][txt3][1])
						print "          polydispersity type: " + txt3 + "; range: " + str(conf_params["model"][txt1][txt2][txt3][0]) + " " + str(conf_params["model"][txt1][txt2][txt3][1])
						
			#now, we create a parameter
			if coupled_name is None:
				if pd_name is None:
					cur_param = param.Param(name=txt1, kind=var_kind, min=min_val, max=max_val, coupled=None, polydispersity=None, pd_min=0.0, pd_max=0.0)
				else:
					cur_param = param.Param(name=txt1, kind=var_kind, min=min_val, max=max_val, coupled=None, polydispersity=txt3, pd_min=pd_minval, pd_max=pd_maxval)
			else:
				if pd_name is None:
					cur_param = param.Param(name=txt1, kind=var_kind, min=min_val, max=max_val, coupled=coupled_name, polydispersity=None, pd_min=0.0, pd_max=0.0)
				else:
					cur_param = param.Param(name=txt1, kind=var_kind, min=min_val, max=max_val, coupled=coupled_name, polydispersity=txt3, pd_min=pd_minval, pd_max=pd_maxval)
			
			#and add it to the array of parameters
			model_params.append(cur_param)
			
		elif params == "Structure_Factor":
			txt1 = str(params)
			for j, sq in enumerate(conf_params["model"][txt1]):
				txt2 = str(sq)
				if txt2 == "name":
					print ""
					print "     Structure factor found:  " + str(conf_params["model"][txt1][txt2])
				else :
					for k, sqpar in enumerate(conf_params["model"][txt1][txt2]):
						txt3 = str(sqpar)
						print "          Structure factor parameter: " + txt2 + "; type: " + txt3 + "; range " + str(conf_params["model"][txt1][txt2][txt3][0]) + " " + str(conf_params["model"][txt1][txt2][txt3][1])
					
	

if __name__=="__main__":
	main()
