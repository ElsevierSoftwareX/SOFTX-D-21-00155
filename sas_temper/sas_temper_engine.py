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

# this handles the set generation
# fconf is the set of configuration parameters for the fitting itself
# modconf is the model to be used and the parameter ranges
# d is the input data to be fit against
def sa_control(fconf, modconf, d):
    # this is the number of refinements to do and how many models to use - yes, they are 'magic numbers'
    refines = 5
    refine_models = 10
    
    # set  up the memory
    res = np.empty(refine_models, "object")            # the parameters returned from the fitting
    mprof = np.empty(refine_models, "object")        # the profiles returned from the fitting
    mprof_usm = np.empty(refine_models, "object")    # the unsmeared profiles returned from the fitting - may be empty
    
    localconf = modelconfig.ModelConfig(modconf.name,modconf.category,modconf.params,modconf.sq)
    
    # st_time = time.time()
    
    # the number of refinement iterations to do
    for j in range(0,refines) :
        for i in range(0,refine_models):
            if j is 0:
                res[i],mprof[i],mprof_usm[i] = sa_engine(fconf,modconf,d)
            else:
                res[i],mprof[i],mprof_usm[i] = sa_engine(fconf,localconf,d)
            
        # refine the ranges to start over
        localconf = sa_refine(refine_models, res)
    
    best = 1000000000.00
    hit = 0
    for k in range(0, refine_models):
        if res[k].chisq < best:
            hit = k
            best = res[k].chisq
    
    #end_time = time.time()
    #dif_time = end_time-st_time
    #junk = "Time for a single model = " + str(dif_time)
    #print(junk)
    
    return res[hit], mprof[hit], mprof_usm[hit]
    

# this is the actual simulated annealing
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
    f = copy.deepcopy(fbest)
    fbest = est_uncerts(d,f,modconf,best_model)

    return fbest, best_model, best_model_usm
    
    
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

# refine the configuration based on a set of results    
def sa_refine(numres, res):
    # make the updated configuration structure and initialize the values
    update_config = modelconfig.ModelConfig(res[0].name,res[0].category,res[0].params,res[0].sq)
    for i, p in enumerate(update_config.params):
        update_config.params[i].min = 1000000000.00
        update_config.params[i].max = -1000000000.00
        
        if p.polydispersity is not None:
            update_config.params[i].polydispersity.min = 1000000000.00
            update_config.params[i].polydispersity.max = -1000000000.00
        
    if update_config.sq is not None:
        for i, sqp in enumerate(update_config.sq.params):
            update_config.sq.params[i].min = 1000000000.00
            update_config.sq.params[i].max = -1000000000.00
            
            if sqp.polydispersity is not None:
                update_config.sq.params[i].polydispersity.min = 1000000000.00
                update_config.sq.params[i].polydispersity.max = -1000000000.00
                
    # first, create the histogram from the current results
    histo_bins = 15
    
    hmin = np.zeros(histo_bins, dtype=np.float64)
    hmax = np.zeros(histo_bins, dtype=np.float64)
    hmid = np.zeros(histo_bins, dtype=np.float64)
    hbin = np.zeros(histo_bins, dtype=np.int)
    chisq = np.zeros(numres, dtype = np.float64)
    
    min = 1000000000.00
    max = -1000000000.00
    for i in range(0,numres):
        chisq[i] = res[i].chisq
        if res[i].chisq < min:
            min = res[i].chisq
        if res[i].chisq > max:
            max = res[i].chisq
    
    if (max/min) < 10.0:
        # use linear bins for the histogram
        min = 0.9*min
        max = 1.1*max
        dx = (max-min)/histo_bins
        
        for i in range(0,histo_bins):
            hmin[i] = i*dx + min
            hmax[i] = (i+1.0)*dx + min
            hmid[i] = (i+0.5)*dx + min
            hbin[i] = 0
            
    else:
        # use logarithmic bins for the histogram
        min = 0.9*m.log10(min)
        max = 1.1*m.log10(max)
        dx = (max-min)/histo_bins
        
        for i in range(0,histo_bins):
            hmin[i] = m.pow(10.0, (i*dx + min))
            hmax[i] = m.pow(10.0, ((i+1.0)*dx + min))
            hmid[i] = m.pow(10.0, ((i+0.5)*dx + min))
            hbin[i] = 0
        
    for i in range(0,numres):
        for j in range(0,histo_bins):
            if (res[i].chisq > hmin[j]) and (res[i].chisq <= hmax[j]):
                hbin[j] += 1
                
    # now that we are done with the histogram, we can 
    # figure out what needs to be updated and how
    # we actually don't need to sort the results
    res_order = np.argsort(chisq)
    check = 0
    hit = 0
    for i in range(0,histo_bins):
        check += hbin[i]
        if hbin[i] > 0:
            hit = i
            break
    check = hbin[hit]
    if check < int(numres/10.0):
        check = int(numres/10.0)
    if check < 3:
        check = 3
        
    # determine the range of values in the set of results
    for i in range(0,check):
        for j, p in enumerate(res[res_order[i]].params):
            if p.val < update_config.params[j].min:
                update_config.params[j].min = p.min
            if p.val > update_config.params[j].max:
                update_config.params[j].max = p.max
                
            if p.polydispersity is not None:
                if p.polydispersity.val < update_config.params[j].polydispersity.min:
                    update_config.params[j].polydispersity.min = p.polydispersity.min
                if p.polydispersity.val > update_config.params[j].polydispersity.max:
                    update_config.params[j].polydispersity.max = p.polydispersity.max
        
        if update_config.sq is not None:
            for j, sqp in enumerate(res[res_order[i]].sq.params):
                if sqp.val < update_config.sq.params[j].min:
                    update_config.sq.params[j].min = sqp.min
                if sqp.val > update_config.sq.params[j].max:
                    update_config.sq.params[j].max = sqp.max
                    
                if sqp.polydispersity is not None:
                    if sqp.polydispersity.val < update_config.sq.params[j].polydispersity.min:
                        update_config.sq.params[j].polydispersity.min = sqp.polydispersity.min
                    if sqp.polydispersity.val > update_config.sq.params[j].polydispersity.max:
                        update_config.sq.params[j].polydispersity.max = sqp.polydispersity.max
        
    return update_config
    
def modify_constraints(analyzed):
    local = modelconfig.ModelConfig(analyzed.name,analyzed.category,analyzed.params,analyzed.sq)
    
    for i, p in enumerate(analyzed.params):
        range = p.max - p.min
        mid = 0.5*(p.min + p.max)
        local.params[i].max = mid + 0.6*range
        local.params[i].min = mid - 0.6*range
        
        if p.polydispersity is not None:
            range = p.polydispersity.max - p.polydispersity.min
            mid = 0.5*(p.polydispersity.min + p.polydispersity.max)
            local.params[i].polydispersity.max = mid + 0.6*range
            local.params[i].polydispersity.min = mid - 0.6*range
        
    if analyzed.sq is not None:
        for i, sqp in enumerate(analyzed.sq.params):
            range = sqp.max - sqp.min
            mid = 0.5*(sqp.min + sqp.max)
            local.sq.params[i].max = mid + 0.6*range
            local.sq.params[i].min = mid - 0.6*range
            
            if sqp.polydispersity is not None:
                range = sqp.polydispersity.max - sqp.polydispersity.min
                mid = 0.5*(sqp.min + sqp.max)
                local.sq.params[i].polydispersity.max = mid + 0.6*range
                local.sq.params[i].polydispersity.min = mid - 0.6*range
    
    return local
    
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
    tprof = sas_data.Model(d, unsmeared = False)
    
    # preparation work for calculating the Jacobian matrix from the derivative
    step = 0.01
    dof = 0
    for i,p in enumerate(modconf.params):
        if p.kind not in ["fixed"]:
            eps.params[i].val = step*(p.max - p.min)
            if eps.params[i].val == 0.0:
                eps.params[i].val = step
            dof = dof + 1
        
        tmp = copy.deepcopy(f)
        tmp.params[i].val = f.params[i].val + eps.params[i].val
        if tmp.params[i].val >= modconf.params[i].max:
            tmp.params[i].val = f.params[i].val - eps.params[i].val
            
        stepped.append(tmp)
        steps.append(eps.params[i].val)
        
    if loc.sq is not None:
        for j, sqp in enumerate(modconf.sq.params):
            if sqp.kind not in ["fixed"]:
                eps.sq.params[j].val = step*(sqp.max-sqp.min)
                if eps.sq.params[j].val == 0.0:
                    eps.sqp.params[j].val = step
                dof = dof + 1
                
            tmp = copy.deepcopy(f)
            tmp.sq.params[j].val = f.sq.params[j].val + eps.sq.params[j].val
            if tmp.sq.params[j].val >= modconf.sq.params[j].max:
                tmp.sq.params[j].val = f.sq.params[j].val - eps.sq.params[j].val
                
            stepped.append(tmp)
            steps.append(eps.sq.params[j].val)
            
    # I don't see a cleaner way to do this...the value of the polydispersity cannot be denoted "fixed"
    for i,p in enumerate(modconf.params):
        if p.polydispersity is not None:
            eps.params[i].polydispersity.val = step*(p.polydispersity.max-p.polydispersity.min)
            if eps.params[i].polydispersity.val == 0.00:
                eps.params[i].polydispersity.val = step
            
            tmp = copy.deepcopy(f)
            tmp.params[i].polydispersity.val = f.params[i].polydispersity.val + eps.params[i].polydispersity.val
            if tmp.params[i].polydispsersity.val >= modconf.params[i].polydispersity.max:
                tmp.params[i].polydispersity.val = f.params[i].polydispersity.val - eps.params[i].polydispersity.val
                
            stepped.append(tmp)
            steps.append(eps.params[i].polydispersity.val)
            
    if loc.sq is not None:
        for j,sqp in enumerate(modconf.sq.params):
            if sqp.polydispersity is not None:
                eps.sq.params[j].polydispersity.val = step*(modconf.sq.params[j].polydispersity.max-modconf.sq.params[j].polydispersity.min)
                if eps.sq.params[j].polydispersity.val == 0.00:
                    eps.sq.params[j].polydispersity.val = step
                    
                tmp = copy.deepcopy(f)
                tmp.sq.params[j].polydispersity.val = f.sq.params[j].polydispersity.val + eps.sq.params[j].polydispersity.val
                if tmp.sq.params[j].polydispersity.val >= modconf.sq.params[j].polydispersity.max:
                    tmp.sq.params[j].polydispersity.val = f.sq.params[j].polydispersity.val - eps.sq.params[j].polydispersity.val
                    
                stepped.append(tmp)
                steps.append(eps.sp.params[j].polydispersity.val)
    
    # and this is where things get ugly
    JT = []
    for w in range(0,len(stepped)):
        #calculate the profiles
        if d.dx is None:
            lprof = sas_calc.calc_profile_usm(d, stepped[w])
        else:
            lprof_usm = sas_calc.calc_profile_usm(d, stepped[w]) 
            lprof = sas_calc.calc_profile(d,stepped[w],lprof_usm)
        
        tprof.y = (lprof.y-best_model.y)/d.dy
        tprof.y = tprof.y*tprof.y
        JT.append(tprof.y/steps[w])
    
    #this is the matrix that we want
    J_T = np.vstack(JT)
    J = J_T.T
    # print(J)
    
    # This is an approximation of the Hessian
    Hess = np.matmul(J_T,J)
    print("Hessian")
    print(Hess)
    
    # Invert it to get the covariance matrix
    Cov = np.linalg.inv(Hess)
    print("Covariance")
    print(Cov)
    
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


