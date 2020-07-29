r"""
sas_temper

polydispersity.py - a class for setting the polydispersity parameters 
                    for a given parameter (i.e. how polydisperse the 
                    radius of a sphere is). After reading in a configuration 
                    file in the YAML format, extracting the information 
                    required to define a model as implemented in the 
                    sasmodels codebase used in SASView and bumps. SASMoodels 
                    requires a name and a set of model-specific parameters to 
                    call the correct code to generate an intensity profile, 
                    which includes specification of polydispersity in some 
                    parameters. The kind of distributions that are acceptable 
                    for sas_temper are those that do not require user input.  
                    The options for a user-specified polydispersity distribution
                    are not acceptable because sas_temper is expected to run
                    in a headless configuration. 

Oak Ridge National Laboratory, 2020

"""

import numpy as np
import math as m

class Polydispersity(object):
	def __init__(self, name, kind, min, max):
		self.name = str(name)		# string:  the name of the parameter
		
		if kind in ["SchulzDispersion","GaussianDispersion","LogNormalDispersion","RectangleDispersion"]:
			self.kind = str(kind)		
		else:
			# throw an error
			err_message = "The polydispersity distribution type " + str(kind) + " is not valid"
			raise Exception(err_message)
		
		self.min = np.float64(min)			# the minimum value
		self.max = np.float64(max)			# the maximum value
		self.val = 0.5*(self.min + self.max)
		self.unc = 0.0
		
		#print ""
		#print "Created a polydispersity parameter " + self.kind + " with a range of widths " + str(self.min) + " to " + str(self.max)
		#print ""
	
