r"""
sas_temper

param.py - model parameters for the simulated annlealing to be
			performed by sas_temper.  The parameters are
			taken from a configuration file provided in the defined
			and parsed using code from this project.  After reading 
			in a configuration file in the format defined as part 
			of the Expert LDRD and extract the information required 
			to define a model as implemented in the sasmodels codebase 
			used in SASView and bumps.  SASMoodels require a name and 
			a set of model-specific parameters to call the correct 
			code to generate an intensity profile. By a parameter, we 
			mean a variable from a model, such as a radius or scattering 
			length density, which fall under a model name given by the 
			"Name" keyword.  The structure factors that exist in SASModels 
			also have parameter sets and the object will be reused there.

Oak Ridge National Laboratory, 2020

"""

import numpy as np
import math as m
import sas_temper.polydispersity as pd

class Param(object):
	def __init__(self, name, kind, min, max, coupled=None, polydispersity=None, pd_min=None, pd_max=None):
		self.name = str(name)		# string:  the name of the parameter
		
		if kind in ["linear","log","integer","coupled","fixed"]:
			self.kind = str(kind)		# string:  linear, log, integer, fixed or coupled
		else:
			# throw an error
			err_message = "The parameter type " + str(kind) + " specified for parameter " + str(name) + " is not valid"
			raise Exception(err_message)
		
		self.min = np.float64(min)			# the minimum value
		self.max = np.float64(max)			# the maximum value
		self.val = 0.5*(self.min + self.max)
		self.unc = 0.0
			
		if coupled is not None:
			self.coupled = str(coupled)
		else:
			self.coupled = None
	
		if polydispersity is not None:
			self.polydispersity = pd.Polydispersity(name=self.name,kind=polydispersity,min=pd_min,max=pd_max)
		else:
			self.polydispersity = None
			
