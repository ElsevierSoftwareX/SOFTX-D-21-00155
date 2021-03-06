r"""
sas_temper

sas_temper_config.py:  code for parsing the YAML config file for sas_temper,
                       which is a YAML 1.2 compliant file containing the
                       information required for the program.  

Oak Ridge National Laboratory, 2020

"""

import yaml

def load_config(config_file="name"):
    with open(config_file, 'r') as f:
        return yaml.load(f, Loader=yaml.SafeLoader)
        

class SAConfiguration(object):
    def __init__(self, params):
        self.datafile = params["name"]
        self.qmin = params["qmin"]
        self.qmax = params["qmax"]
        self.output = params["output"]
        self.temperatures = params["temperatures"]
        self.temp_rate = params["temperature_rate"]
        self.param_rate = params["parameter_rate"]
        self.iterations = params["iterations"]
        self.models = params["models"]
        
        #print ("Set up the SA configuration parameter class")
        #print ("")
        #print ("The file is " + str(self.datafile))
        #outbuf = "    to be fit over q = " + str(self.qmin) + " to " + str(self.qmax)
        #print(outbuf)
        #outbuf = "The output filename base is " + str(self.output)
        #print(outbuf)
        #print ("")
        #print ("The parameters for the simulated annealing:")
        #outbuf = "     temperatures: " + str(self.temperatures)
        #print(outbuf)
        #outbuf = "     temperature_rate: " + str(self.temp_rate)
        #print(outbuf)
        #outbuf = "     parameter_rate: " + str(self.param_rate)
        #print(outbuf)
        #outbuf = "     iterations: " + str(self.iterations)
        #print(outbuf)
        #outbuf = "     models to generate: " + str(self.models)
        #print(outbuf)
        
