r"""
sas_temper

modelconfig.py - a class that contains information about the model
                to be used.  It is initiated by reading the config
                file.  This container class is also used during the 
                simulated annealing for holding the parameters and
                ranges and the like.  it also includes a function
                for dealing with log and coupled parameters for 
                the profile calculation and output. Structure factors
                are parts of models, even though SASModels treats them as models.

Oak Ridge National Laboratory, 2020

"""

import math as m
import copy

import sas_temper.param as param
import sas_temper.structurefactor as structurefactor

class ModelConfig(object):
    def __init__(self, name, category, parameters, str_factor=None):
        self.name = name
        self.category = category
        self.chisq = 1000000000.00   # determined from the values contained in the 'val of each parameter
        
        self.params = []
        self.free_params = 0
        for a, par in enumerate(parameters):
            self.params.append(par)
            self.free_params += 1
            
        if str_factor is not None:
            self.sq = structurefactor.StructureFactor(str_factor.type,str_factor.params)
            for b, spq in enumerate(str_factor.params):
                self.free_params += 1
        else:
            self.sq = None
        

# convert a config with log and coupled parameters to one with real values
# for use in the calculations and in the outputting of results.
def convert_conf(mc):
    local = ModelConfig(mc.name, mc.category, mc.params, mc.sq)
    
    # need to set this for using the code in output.py
    local = copy.deepcopy(mc)
    
    # not pretty, but it has to be done
    for i, p in enumerate(local.params):
        # first, we take care of the log values
        if p.kind=="log":
            local.params[i].val = m.pow(10.0,p.val)
            #local.params[i].unc = 0.5*(m.pow(10.0,(p.val+p.unc))-m.pow(10.0,(p.val-p.unc)))
        
        # and do the coupling while we're here - just ugly, but whatever
        if p.coupled is not None:
            for j,q in enumerate(local.params):
                if q.name==p.coupled:
                    local.params[i].val = p.val*q.val
                    #local.params[i].unc = p.val*q.unc
            if local.sq is not None:
                for k,r in enumerate(local.sq.params):
                    if r.name==p.coupled:
                        local.params[i].val = p.val*r.val
                        #local.params[i].unc = p.val*r.unc
                    
    if local.sq is not None:
        for i,sqp in enumerate(local.sq.params):
            # first, we take care of the log values
            if sqp.kind=="log":
                local.sq.params[i].val = m.pow(10.0,sqp.val)
                #local.sq.params[i].unc = 0.5*(m.pow(10.0,(sqp.val+sqp.unc))-m.pow(10.0,(sqp.val-sqp.unc)))
            
            # and do the coupling while we're here - just ugly, but whatever
            if sqp.coupled is not None:
                for j,q in enumerate(local.params):
                    if q.name==sqp.coupled:
                        local.sq.params[i].val = sqp.val*q.val
                        #local.sq.params[i].unc = sqp.val*q.unc
                if local.sq is not None:
                    for k,r in enumerate(local.sq.params):
                        if r.name==sqp.coupled:
                            local.sq.params[i].val = sqp.val*r.val
                            #local.sq.params[i].val = sqp.val*r.unc
    
    return local
