r"""
sas_temper

parse_conf.py:  code for parsing the YAML configuration file.  
                Read in a YAML 1.2 compliant configuration file
                for the fitting engine that is built on simulated 
                annealing.  

Oak Ridge National Laboratory, 2020

"""

import sys
import sas_temper.modelconfig as modelconfig
import sas_temper.sas_temper_config as config
import sas_temper.param as param
import sas_temper.polydispersity as polydispersity
import sas_temper.structurefactor as structurefactor

def parse_config(name):
	conf_params = config.load_config(config_file=name)
	
	sa_parameters = {}
	sa_parameters["name"] = str(conf_params["file"]["name"])
	sa_parameters["qmin"] = float(conf_params["file"]["qmin"])
	sa_parameters["qmax"] = float(conf_params["file"]["qmax"])
	sa_parameters["output"] = str(conf_params["output_files"])
	sa_parameters["temperatures"] = float(conf_params["sa_parameters"]["temperatures"])
	sa_parameters["temperature_rate"] = float(conf_params["sa_parameters"]["temperature_rate"])
	sa_parameters["parameter_rate"] = float(conf_params["sa_parameters"]["parameter_rate"])
	sa_parameters["iterations"] = int(conf_params["sa_parameters"]["iterations"])
	sa_parameters["models"] = int(conf_params["sa_parameters"]["models_to_generate"])
	
	sasTemperConf = config.SAConfiguration(sa_parameters)
	
    model_name = None
    model_category = None
	model_params = []
	s_of_q = None
	s_of_q_params = []
	
	for i, params in enumerate(conf_params["model"]):
		if params not in ["name","category","Category","Structure_Factor"]:
			txt1 = str(params)
			#print ""
			#print "     parameter name:  " + txt1
			
			coupled_name=None
			pd_name=None
			var_kind=None
			min_val = 1.0
			max_val = 1.0
			pd_minval = 1.0
			pd_maxval = 1.0
			for j, mod in enumerate(conf_params["model"][txt1]):
				txt2 = str(mod)
				if mod in ["fixed", "linear", "log", "integer"]:
					var_kind = txt2
					if txt2 in ["linear", "log", "integer"]:
						#print "          type: " + txt2 + "; value: " + str(conf_params["model"][txt1][txt2][0])
						min_val = float(conf_params["model"][txt1][txt2][0])
						max_val = float(conf_params["model"][txt1][txt2][1])
					else:
						#print "          type: " + txt2 + "; range: " + str(conf_params["model"][txt1][txt2][0]) + " " + str(conf_params["model"][txt1][txt2][1])
						min_val = float(conf_params["model"][txt1][txt2][0])
						max_val = float(conf_params["model"][txt1][txt2][0])
				elif mod == "coupled":
					coupled_name = str(conf_params["model"][txt1][txt2][0])
					#print "          coupled to: " + str(conf_params["model"][txt1][txt2][0])
				elif mod == "polydispersity":
					for k, poly in enumerate(conf_params["model"][txt1][txt2]):
						txt3 = str(poly)
						pd_name = txt3
						pd_minval = float(conf_params["model"][txt1][txt2][txt3][0])
						pd_maxval = float (conf_params["model"][txt1][txt2][txt3][1])
						#print "          polydispersity type: " + txt3 + "; range: " + str(conf_params["model"][txt1][txt2][txt3][0]) + " " + str(conf_params["model"][txt1][txt2][txt3][1])
						
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
			
		elif params in ["Structure_Factor"]:
			txt1 = str(params)
			for j, sq in enumerate(conf_params["model"][txt1]):
				txt2 = str(sq)
				sq_coupled_name=None
				sq_pd_name=None
				sq_var_kind=None
				sq_min_val = 1.0
				sq_max_val = 1.0
				sq_pd_min_val = 1.0
				sq_pd_max_val = 1.0
				if txt2 in ["name"]:
					sq_name = str(conf_params["model"][txt1][txt2])
				else :
					sq_coupled_name=None
					sq_pd_name=None
					sq_var_kind=None
					for j, sq_mod in enumerate(conf_params["model"][txt1][txt2]):
						txt3 = str(sq_mod)
						if sq_mod in ["fixed", "linear", "log", "integer"]:
							sq_var_kind = txt3
							if txt3 in ["linear", "log", "integer"]:
								sq_min_val = float(conf_params["model"][txt1][txt2][txt3][0])
								sq_max_val = float(conf_params["model"][txt1][txt2][txt3][1])
							else:
								sq_min_val = float(conf_params["model"][txt1][txt2][txt3][0])
								sq_max_val = float(conf_params["model"][txt1][txt2][txt3][0])
						elif sq_mod in ["coupled"]:
							sq_coupled_name = str(conf_params["model"][txt1][txt2][txt3][0])
						elif sq_mod in ["polydispersity"]:
							for k, sq_poly in enumerate(conf_params["model"][txt1][txt2][txt3]):
								txt4 = str(sq_poly)
								sq_pd_name = txt4
								sq_pd_minval = float(conf_params["model"][txt1][txt2][txt3][txt4][0])
								sq_pd_maxval = float(conf_params["model"][txt1][txt2][txt3][txt4][1])
						
					#now, we create a parameter
					if sq_coupled_name is None:
						if sq_pd_name is None:
							sq_cur_param = param.Param(name=txt2, kind=sq_var_kind, min=sq_min_val, max=sq_max_val, coupled=None, polydispersity=None, pd_min=0.0, pd_max=0.0)
						else:
							sq_cur_param = param.Param(name=txt2, kind=sq_var_kind, min=sq_min_val, max=sq_max_val, coupled=None, polydispersity=txt3, pd_min=pd_minval, pd_max=pd_maxval)
					else:
						if sq_pd_name is None:
							sq_cur_param = param.Param(name=txt2, kind=sq_var_kind, min=sq_min_val, max=sq_max_val, coupled=sq_coupled_name, polydispersity=None, pd_min=0.0, pd_max=0.0)
						else:
							sq_cur_param = param.Param(name=txt2, kind=sq_var_kind, min=sq_min_val, max=sq_max_val, coupled=sq_coupled_name, polydispersity=txt4, pd_min=sq_pd_minval, pd_max=sq_pd_maxval)
	
					#and add it to the array of structure factor parameters
					s_of_q_params.append(sq_cur_param)
				
			#we need to actually create a structure factor before we make the final model configuration	
			s_of_q = structurefactor.StructureFactor(sq_name,s_of_q_params)
			
		elif params in ["name"]:
			model_name = str(conf_params["model"]["name"])
		elif params in ["Category"]:
			model_category = str(conf_params["model"]["Category"])
	
	modelconf = modelconfig.ModelConfig(model_name,model_category,model_params,str_factor=s_of_q)
	
	return sasTemperConf, modelconf
	
