r"""
sas_temper

structurefactor.py - a class that specifies a structure factor for a 
                     model.  A YAML configuration file for sas_temper
                     can also be used to specify that a model be multiplied
                     by a structure factor.  In a sasmodels, a structure 
                     factor is just another form of model for the purpose of
                     the calculation.  Once calculated, sas_temper multiplies
                     F(q)*S(q) for I(q) to compare against the data.  

Oak Ridge National Laboratory, 2020

"""

import sas_temper.param as param

class StructureFactor(object):
    def __init__(self, type, params):
        self.params = []
        
        if type in ["harsphere", "hayter_msa", "squarewell", "stickyhardsphere"]:
            self.type = type
        else:
            # we just set this rather then throw an error
            self.type = "hardsphere"
        
        for a, par in enumerate(params):
            self.params.append(par)
            
        
    
