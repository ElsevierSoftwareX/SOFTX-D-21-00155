r"""
sas_temper

sas_temper_engine.py:  The actual simulated annealing algorithm is 
                       here, as is the code for refinement process 
                       and the wrapper function that controls both.

Oak Ridge National Laboratory, 2020

"""

import numpy as np
import math as m
import copy
import time

import sas_temper.modelconfig as modelconfig
import sas_temper.sas_temper_config as config
import sas_temper.sas_data as sas_data
import sas_temper.sas_calc as sas_calc

# this function handles the set generation
# fconf is the set of configuration parameters for the fitting itself
# modconf is the model to be used and the parameter ranges
# d is the input data to be fit against
def sa_control(fconf, modconf, d):
    # we can get right down to business since this is a single minimization.
    
    st_time = time.time()
    
    res,mprof,mprof_usm = sa_engine(fconf,modconf,d)
    
    end_time = time.time()
    dif_time = end_time-st_time
    junk = "Time for a single model = " + str(dif_time)
    print(junk)
    
    return res, mprof, mprof_usm
    

# this function performs the actual simulated annealing
# fconf is the set of configuration parameters for the fitting itself
# modconf is the model to be used and the parameter ranges
# d is the input data to be fit against
def sa_engine(fconf, modconf, d): 
    # set up the results containers
    fbest = modelconfig.ModelConfig(modconf.name,modconf.category,modconf.params,modconf.sq)
    cur = modelconfig.ModelConfig(modconf.name,modconf.category,modconf.params,modconf.sq)
    f = modelconfig.ModelConfig(modconf.name,modconf.category,modconf.params,modconf.sq)
    
    # these are the model profiles and the model
    # profiles to keep for the final output
    model_usm = sas_data.Model(d, unsmeared = True)
    best_model_usm = sas_data.Model(d, unsmeared = True)
    model = sas_data.Model(d, unsmeared = False)
    best_model = sas_data.Model(d, unsmeared = False)
    
    # initialize things outside the while loop
    f = define_model(1,modconf,10.0,1.0,cur)
    fbest = copy.deepcopy(f)
    cur = copy.deepcopy(f)
    #calculate the profiles
    if d.dx is None:
        model = sas_calc.calc_profile_usm(d,f)
    else:
        model_usm = sas_calc.calc_profile_usm(d,f) 
        model = sas_calc.calc_profile(d,f,model_usm)
    
    f.chisq = sas_calc.chisq(f,d,model)
    fbest.chisq = f.chisq
    cur.chisq = f.chisq
    bchi = f.chisq
    
    # this is the main simulated annealing loop
    temp = 10.0
    r = 1.0
    schedule = 1
    hit = False
    while schedule <= fconf.temperatures:
        iters = 1
        
        while iters <= fconf.iterations:
            #define the model to test
            f = define_model(schedule,modconf,temp,r,cur)
            
            #calculate the profiles
            if d.dx is None:
                model = sas_calc.calc_profile_usm(d,f)
            else:
                model_usm = sas_calc.calc_profile_usm(d,f)
                model = sas_calc.calc_profile(d,f,model_usm)
                
            f.chisq = sas_calc.chisq(f,d,model)
            
            if f.chisq < fbest.chisq:
                hit = True
            
                # we have found our new best overall
                if f.chisq < bchi:
                    bchi = f.chisq
                    
                    fbest = copy.deepcopy(f)
                    best_model = copy.deepcopy(model)
                    best_model_usm = copy.deepcopy(model_usm)
                
            else:
                energy = f.chisq - cur.chisq
                val = m.exp(-1.0*energy/temp)
                if frand(0.0,1.0) < val:
                    hit = True
                else:
                    hit = False
                
            if hit:
                cur = copy.deepcopy(f)
                
            iters = iters + 1
            
            #noise = "schedule " + str(schedule) + "; iteration " + str(iters)
            #print(noise)
        
        # we drop the temperature, tighten the range and increase schedule
        temp = temp*fconf.temp_rate
        r = r*fconf.param_rate
        schedule = schedule + 1
        
    
    # as the final step, we estimate the uncertainties with the jacobian
    fbest_unc = est_uncerts(d,fbest,modconf,best_model)

    return fbest_unc, best_model, best_model_usm
    
    
# the python version of my tried and true c code
def frand(min, max):
    val = (max-min)*np.random.random() + min
    
    return val


# define a new random model, possibly from the current state
def define_model(schedule, modconf, temperature, rval, current):
    local = modelconfig.ModelConfig(modconf.name,modconf.category,modconf.params,modconf.sq)
    
    if schedule is 1:
        for i, p in enumerate(modconf.params):
            if p.kind in ["integer"]:
                local.params[i].val = int(frand(p.min,p.max+1.0))
            else:
                local.params[i].val = frand(p.min,p.max)
            
            if p.polydispersity is not None:
                # polydispersity is a funny sort of parameter and is never anything but linear
                local.params[i].polydispersity.val = frand(p.polydispersity.min,p.polydispersity.max)
            
        if local.sq is not None:
            for i, sqp in enumerate(modconf.sq.params):
                if sqp.kind in ["integer"]:
                    local.sq.params[i].val = int(frand(sqp.min,sqp.max+1.0))
                else:
                    local.sq.params[i].val = frand(sqp.min,sqp.max)
                
                if sqp.polydispersity is not None:
                    # polydispersity is a funny sort of parameter and is never anything but linear
                    local.sq.params[i].polydispersity.val = frand(sqp.polydispersity.min,sqp.polydispersity.max)
    else :
        # pick a value from within the constrained range
        for i, p in enumerate(current.params):
            move = frand(-0.1*rval*p.val,0.1*rval*p.val)
            if p.kind in ["integer"]:
                local.params[i].val = int(p.val + move)
            else:
                local.params[i].val = p.val + move
                
            if local.params[i].val>modconf.params[i].max:
                local.params[i].val = modconf.params[i].max
            if local.params[i].val<=modconf.params[i].min:
                local.params[i].val = modconf.params[i].min
            
            if p.polydispersity is not None:
                move = frand(-0.1*rval*p.polydispersity.val,0.1*rval*p.polydispersity.val)
                #catch the integer types
                if p.kind in ["integer"]:
                    local.params[i].polydispersity.val = int(p.polydispersity.val + move)
                else:
                    local.params[i].polydispersity.val = p.polydispersity.val + move
                
                if local.params[i].polydispersity.val > modconf.params[i].polydispersity.max:
                    local.params[i].polydispersity.val = modconf.params[i].polydispersity.max
                if local.params[i].polydispersity.val < modconf.params[i].polydispersity.min:
                    local.params[i].polydispersity.val = modconf.params[i].polydispersity.min    
        
        if local.sq is not None:        
            for i, sqp in enumerate(current.sq.params):
                move = frand(-0.1*rval*sqp.val,0.1*rval*sqp.val)
                if sqp.kind in ["integer"]:
                    local.sq.params[i].val = int(sqp.val + move)
                else:
                    local.sq.params[i].val = sqp.val + move
                
                if local.sq.params[i].val>modconf.sq.params[i].max:
                    local.sq.params[i].val = modconf.sq.params[i].max
                if local.sq.params[i].val<=modconf.sq.params[i].min:
                    local.sq.params[i].val = modconf.sq.params[i].min
                    
                if sqp.polydispersity is not None:
                    move = frand(-0.1*rval*sqp.polydispersity.val,0.1*rval*sqp.polydispersity.val)
                    if sqp.kind in ["integer"]:
                        local.sq.params[i].polydispersity.val = int(sqp.polydispersity.val + move)
                    else:
                        local.sq.params[i].polydispersity.val = sqp.polydispersity.val + move
                        
                    if local.sq.params[i].polydispersity.val > modconf.sq.params[i].polydispersity.max:
                        local.sq.params[i].polydispersity.val = modconf.sq.params[i].polydispersity.max
                    if local.sq.params[i].polydispersity.val < modconf.sq.params[i].polydispersity.min:
                        local.sq.params[i].polydispersity.val = modconf.sq.params[i].polydispersity.min
    
    return local

# perform the traditional estimation of the uncertainties in the
# fitting parameters by looking at the partial derivatives using the
# approach that is used in Paul Kienzle's "bump", as is noted below.    
def est_uncerts(d, f, modconf, best_model):
    # this is our ugly way of getting at this matrix of derivatives
    loc = modelconfig.ModelConfig(f.name,f.category,f.params,f.sq)
    loc = copy.deepcopy(f)
    eps = modelconfig.ModelConfig(f.name,f.category,f.params,f.sq)
    tmp = modelconfig.ModelConfig(f.name,f.category,f.params,f.sq)
    stepped = []
    steps = []
    
    # local profile for the calculation of the derivative
    lprof_usm = sas_data.Model(d, unsmeared = True)
    lprof = sas_data.Model(d, unsmeared = False)
    
    # preparation work for calculating the Jacobian matrix from the derivative
    step = 0.001
    for i,p in enumerate(modconf.params):
        if p.kind not in ["fixed"]:
            eps.params[i].val = step*(p.max - p.min)
            if eps.params[i].val == 0.0:
                eps.params[i].val = step
        else: 
            eps.params[i].val = 0.01*step
        
        tmp = copy.deepcopy(loc)
        tmp.params[i].val = loc.params[i].val + eps.params[i].val
        # print("tmp.params[i].val = "+str(tmp.params[i].val)+"   f.params[i].val = "+str(loc.params[i].val)+"   eps.params[i].val = "+str(eps.params[i].val))
        if tmp.params[i].val >= modconf.params[i].max:
            tmp.params[i].val = loc.params[i].val - eps.params[i].val
            
        stepped.append(tmp)
        steps.append(eps.params[i].val)
        
    if loc.sq is not None:
        for j, sqp in enumerate(modconf.sq.params):
            if sqp.kind not in ["fixed"]:
                eps.sq.params[j].val = step*(sqp.max-sqp.min)
                if eps.sq.params[j].val == 0.0:
                    eps.sqp.params[j].val = step
            else: 
                eps.sqp.params[j].val = 0.01*step
                
            tmp = copy.deepcopy(loc)
            tmp.sq.params[j].val = loc.sq.params[j].val + eps.sq.params[j].val
            if tmp.sq.params[j].val >= modconf.sq.params[j].max:
                tmp.sq.params[j].val = loc.sq.params[j].val - eps.sq.params[j].val
                
            stepped.append(tmp)
            steps.append(eps.sq.params[j].val)
            
    # I don't see a cleaner way to do this...the value of the polydispersity cannot be denoted "fixed"
    for i,p in enumerate(modconf.params):
        if p.polydispersity is not None:
            eps.params[i].polydispersity.val = step*(p.polydispersity.max-p.polydispersity.min)
            if eps.params[i].polydispersity.val == 0.00:
                eps.params[i].polydispersity.val = step
            
            tmp = copy.deepcopy(loc)
            tmp.params[i].polydispersity.val = loc.params[i].polydispersity.val + eps.params[i].polydispersity.val
            if tmp.params[i].polydispersity.val >= modconf.params[i].polydispersity.max:
                tmp.params[i].polydispersity.val = loc.params[i].polydispersity.val - eps.params[i].polydispersity.val
                
            stepped.append(tmp)
            steps.append(eps.params[i].polydispersity.val)
            
    if loc.sq is not None:
        for j,sqp in enumerate(modconf.sq.params):
            if sqp.polydispersity is not None:
                eps.sq.params[j].polydispersity.val = step*(modconf.sq.params[j].polydispersity.max-modconf.sq.params[j].polydispersity.min)
                if eps.sq.params[j].polydispersity.val == 0.00:
                    eps.sq.params[j].polydispersity.val = step
                    
                tmp = copy.deepcopy(loc)
                tmp.sq.params[j].polydispersity.val = loc.sq.params[j].polydispersity.val + eps.sq.params[j].polydispersity.val
                if tmp.sq.params[j].polydispersity.val >= modconf.sq.params[j].polydispersity.max:
                    tmp.sq.params[j].polydispersity.val = loc.sq.params[j].polydispersity.val - eps.sq.params[j].polydispersity.val
                    
                stepped.append(tmp)
                steps.append(eps.sp.params[j].polydispersity.val)
    
    # and this is where things get ugly
    JT = []
    for w in range(0,len(stepped)):
        #for a,m in enumerate(stepped[w].params):
        #    print("stepped["+str(w)+"] parameters " + str(m.val))
            
        #calculate the profiles
        if d.dx is None:
            lprof = sas_calc.calc_profile_usm(d, stepped[w])
        else:
            lprof_usm = sas_calc.calc_profile_usm(d, stepped[w]) 
            lprof = sas_calc.calc_profile(d,stepped[w],lprof_usm)
        
        tprof = sas_data.Model(d, unsmeared = False)
        for z in range(0,len(tprof.y)):
            tprof.y[z] = 0.5*(lprof.y[z]-best_model.y[z])/(d.dy[z]*steps[w])
            #if z==0:
                #print(str(best_model.y[z]) + "     " + str(lprof.y[z]) + "     " + str(steps[w]) + "    tprof.y[0] = "+str(tprof.y[z]))
            
        JT.append(tprof.y)
    
    #this is the matrix that we want
    J_T = np.vstack(JT)
    J = J_T.T
    # print("Jacobian")
    # print(J)
    
    # this is the singlular value decomposition of J
    # see https://github.com/bumps/bumps/blob/master/bumps/lsqerror.py lines 237-239
    # The tolerance shown cuts down on singlular values
    # bumps as a whole may not use this directly.  have to give it a try, though
    U, S, V = np.linalg.svd(J, 0)
    tol = 1e-8
    S[S <= tol] = tol
    
    # This is the work around to avoid the direct inversion of the psuedo-Hessian
    Cov = np.dot(V.T.conj() / S**2, V)
    
    # the diagonal should be only as long as the number of parameters
    errs = np.sqrt(np.diag(Cov))

    # these final values need to be inserted into the structure to be returned
    # note that fixed parameters have their uncertainty set to 0.0 here to avoid
    # the use of the tolerance applied above that avoids dividing by zero
    k = 0
    for i,p in enumerate(loc.params):
        if loc.params[i].kind not in ["fixed"]:
            loc.params[i].unc = errs[k]
        else:
            loc.params[i].unc = 0.00
        k = k + 1
    if loc.sq is not None:
        for j, sqp in enumerate(loc.sq.params):
            if loc.sq.params[j].kind not in ["fixed"]:
                loc.sq.params[j].unc = errs[k]
            else:
                loc.sq.params[j].unc = 0.00
            k = k + 1
    # we filled it like this above...
    for i,p in enumerate(loc.params):
        if p.polydispersity is not None:
            loc.params[i].polydispersity.unc = errs[k]
            k = k + 1
    if loc.sq is not None:
        for j,sqp in enumerate(loc.sq.params):
            if sqp.polydispersity is not None:
                loc.sq.params[j].polydispersity.unc = errs[k]
                k = k + 1
    
    # finally, we can return from this
    return(loc)


