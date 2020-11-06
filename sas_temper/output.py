r"""
sas_temper

output.py:  code for writing fitting results and the results of
            the analysis performed on the set of results. This code 
            does two different things.  It outputs the model intensity 
            profiles with parameters found during the fitting.  It 
            also outputs the final table of averages, standard deviations 
            and the summary table of the models found during the fitting.

Oak Ridge National Laboratory, 2020

"""

from matplotlib import pyplot as plt
import numpy as np

import sas_temper.sas_temper_config as config
import sas_temper.sas_data as sas_data
import sas_temper.modelconfig as modelconfig
import sas_temper.param as param
import sas_temper.structurefactor as structurefactor
import sas_temper.polydispersity as polydispersity

#this function outputs a single model profile with the associated fit parameters and uncertainties
def outputSingleRes(conf, d, m, mnum, res):
    # convert the config to real units
    lres = modelconfig.convert_conf(res)

    name = str(conf.output)+"%02d.txt" %(mnum)
    with open(name, 'w') as f:
        # a little header information
        buff = "# input data: " + str(conf.datafile) + "\n"
        f.write(buff)
        buff = "# the name of the model used:  " + str(lres.name) + "\n"
        f.write(buff)
        buff = "# chi-squared is %6.6f \n" %(lres.chisq)
        f.write(buff)
        buff = "# list of model fitting parameters and values\n"
        f.write(buff) 
        
        # these are the parameters and values
        for i, p in enumerate(lres.params):
            if p.polydispersity is None:
                if p.kind not in ["fixed"]:
                    if p.kind in ["integer"]:
                        buff = "# "+str(p.name)+ " = %6d +/- %6.6f\n" %(p.val, p.unc)
                    else:
                        buff = "# "+str(p.name)+ " = %6.6f +/- %6.6f\n" %(p.val, p.unc)
                else:
                    buff = "# "+str(p.name)+ " = %6.6f; fixed\n" %(p.val)
            else: 
                if p.kind not in ["fixed"]:
                    if p.kind in ["integer"]:
                        buff = "# "+str(p.name)+ " = %6d +/- %6.6f, polydispersity %s = %6.6f +/- %6.6f\n" %(p.val, p.unc, str(p.polydispersity.kind), p.polydispersity.val, p.polydispersity.unc)
                    else:
                        buff = "# "+str(p.name)+ " = %6.6f +/- %6.6f, polydispersity %s = %6.6f +/- %6.6f\n" %(p.val, p.unc, str(p.polydispersity.kind), p.polydispersity.val, p.polydispersity.unc)
                else:
                    # the polydispersity parameter cannot be "fixed" at present
                    buff = "# "+str(p.name)+ " = %6.6f; fixed, polydispersity %s = %6.6f +/- %6.6f\n" %(p.val, str(p.polydispersity.kind), p.polydispersity.val, p.polydispersity.unc)
            f.write(buff)
            
            # in case this was coupled
            if p.coupled is not None:
                buff = "# \tcoupled to " + str(p.coupled) + "\n"
                f.write(buff)
        
        # The structure factor applied is also a parameter and is included here        
        if lres.sq is not None:
            buff = "# Structure factor applied: " + str(lres.sq.type) + "\n"
            f.write(buff)
            for i, sqp in enumerate(lres.sq.params):
                if sqp.polydispersity is None:
                    if sqp.kind not in ["fixed"]:
                        if sqp.kind in ["integer"]:
                            buff = "# " + str(sqp.name) + " = %6d +/- %6.6f \n" %(sqp.val, sqp.unc)
                        else:
                            buff = "# " + str(sqp.name) + " = %6.6f +/- %6.6f \n" %(sqp.val, sqp.unc)
                    else:
                        buff = "# " + str(sqp.name) + " = %6.6f; fixed \n" %(sqp.val)
                else: 
                    if sqp.kind not in ["fixed"]:
                        if sqp.kind in ["integer"]:
                            buff = "# " + str(sqp.name) + " = %6d +/- %6.6f, polydispersity %s = %6.6f +/- %6.6f\n" %(sqp.val, sqp.unc, str(sqp.polydispersity.kind), sqp.polydispersity.val, sqp.polydispersity.unc)
                        else:
                            buff = "# " + str(sqp.name) + " = %6.6f +/- %6.6f, polydispersity %s = %6.6f +/- %6.6f\n" %(sqp.val, sqp.unc, str(sqp.polydispersity.kind), sqp.polydispersity.val, sqp.polydispersity.unc)
                    else:
                        # the polydispersity parameter cannot be "fixed" at present
                        buff = "# " + str(sqp.name) + " = %6.6f; fixed, polydispersity %s = %6.6f +/- %6.6f\n" %(sqp.val, str(sqp.polydispersity.kind), sqp.polydispersity.val, sqp.polydispersity.unc)
                f.write(buff)
                
                # in case we have a coupled parameter
                if sqp.coupled is not None :
                    buff = "# \tcoupled to " + str(sqp.coupled) + "\n"
                    f.write(buff)
                
        
        buff = "# the model profile q and I(q)\n"
        for i in range(0,len(m.x)):
            buff = "%6.6f\t%6.6f\n" %(m.x[i], m.y[i]) 
            f.write(buff)
    
    # end by closing the file
    f.close()
    
    # output a plot of the data and the fit curve
    outputFitCurve(conf,d,m,mnum,lres.chisq)

    
# this function takes the set of results, calculates the average and standard deviations, and outputs them
def outputSetRes(conf, res):
    #determine how many things we will be doing statistics on
    parms = 0
    for i, p in enumerate(res[0].params):
        if p.polydispersity is None:
            parms = parms + 1
        else: 
            parms = parms + 2
    
    if res[0].sq is not None:
        for i, sqp in enumerate(res[0].sq.params):
            if sqp.polydispersity is None:
                parms = parms + 1
            else: 
                parms = parms + 2
    
    lres = np.empty(conf.models, 'object')
    chisq = np.empty(conf.models,dtype='float64')
    names = np.empty(parms,dtype='object')
    varkinds = np.empty(parms,dtype='object')
    vals = np.empty((parms,conf.models),dtype='float64')
    for w in range(0,conf.models):
        # convert the results to real units
        lres[w] = modelconfig.convert_conf(res[w])
        
        chisq[w] = res[w].chisq
        j = 0
        while j<parms:
            for i, p in enumerate(lres[w].params):
                if p.polydispersity is None:
                    names[j] = p.name
                    varkinds[j] = p.kind
                    vals[j][w] = p.val
                    j+=1
                else: 
                    names[j] = p.name
                    varkinds[j] = p.kind
                    vals[j][w] = p.val
                    j = j + 1
                    names[j] = p.name+".pd"
                    vals[j][w] = p.polydispersity.val
                    j = j + 1
        
            if lres[w].sq is not None:
                for i, sqp in enumerate(lres[w].sq.params):
                    if sqp.polydispersity is None:
                        names[j] = sqp.name
                        varkinds[j] = sqp.kind
                        vals[j][w] = sqp.val
                        j = j + 1
                    else: 
                        names[j] = sqp.name
                        varkinds[j] = sqp.kind
                        vals[j][w] = sqp.val
                        j = j + 1
                        names[j] = sqp.name+".pd"
                        vals[j][w] = sqp.polydispersity.val
                        j = j + 1

    #calculate the averages and standard deviations of the values
    ave = np.empty(parms,dtype=np.float64)
    std = np.empty(parms,dtype=np.float64)
    cor = np.empty([parms,parms],dtype=np.float64)
    for i in range(0,parms):
        ave[i] = np.average(vals[i])
        std[i] = np.std(vals[i])
    
    #this is to address a bad division that does not crash the code
    with np.errstate(divide='ignore', invalid='ignore'):
        cor = np.corrcoef(vals,y=None,rowvar=True)
        
    #this gets rid of junk values when parameters are fixed
    for i in range(0,parms):
        for j in range(0,parms):
            if varkinds[i] in ["fixed"]:
                cor[i][j] = 1.0
            if varkinds[j] in ["fixed"]:
                cor[i][j] = 1.0
    
    #now, we can write the analysis of the results and the list of results
    name = str(conf.output)+"_set_analysis.txt" 
    with open(name, 'w') as f:
        buff = "# analysis of the set of models " + conf.output + "\n"
        f.write(buff)
        for j in range(0,parms):
            if varkinds[j] in ["fixed"]:
                buff = "# " + names[j] + " = %6.6f; fixed\n" %(ave[j])
            else:
                buff = "# " + names[j] + " = %6.6f +/- %6.6f\n" %(ave[j], std[j])
            f.write(buff)
            
        buff = "# Pearson product-moment correlation coefficients\n"
        f.write(buff)
        for i in range(0,parms):
            buff = "# "
            for j in range(0,parms):
                buff += "%6.6f\t" %(cor[i][j])
            buff += "\n" 
            f.write(buff)
            
        buff = "# chi-squared and parameters from the individual models, in the order of the list of averages above\n"
        f.write(buff)
        for i in range(0,conf.models):
            buff = ""
            buff += "%6.6f\t" %(chisq[i])
            for j in range(0,parms):
                if varkinds[j] in ["integer"]:
                    buff += "%6d\t" %(vals[j][i])
                else:
                    buff += "%6.6f\t" %(vals[j][i])
            buff += "\n"
            f.write(buff)
            
    f.close()
    
    # this outputs a set of graphs showing the relationships between all pairs of variables
    # these are simple scatter plots of one variable plotted against the other, but it should be useful
    for i in range(0,(len(names)-1)):
        if varkinds[i] not in ["fixed"]:
            for j in range((i+1),len(names)):
                if varkinds[j] not in ["fixed"]:
                    fig = plt.figure()
                    grph = fig.add_subplot(1,1,1)
                    grph.set_autoscale_on(True)
                    grph.plot(vals[i],vals[j],'ko')
                    grph.set_title(str(names[j])+" vs. "+str(names[i]))
                    
                    oname = str(conf.output)+"_"+str(names[j])+"_"+str(names[i])+".png"
                    fig.savefig(oname,format="png")
                    plt.close(fig)
                    
    # output histograms of the non-fixed parameters
    for i in range(0,len(names)):
        if varkinds[i] not in ["fixed"]:
            fig = plt.figure(figsize = [4,4], dpi=100)
            grph = fig.add_subplot(1,1,1)
            grph.set_autoscale_on(True)
            grph.hist(vals[i],bins=5,color = 'r', rwidth=0.9)
            grph.set_title("Histogram of "+str(names[i]))
            
            oname = str(conf.output)+"_"+str(names[i])+"_histogram.png"
            fig.savefig(oname,format="png")
            plt.close(fig)
            
    # output plots of the values of the parameters vs. chi-squared
    for i in range(0, len(names)):
        if varkinds[i] not in ["fixed"]:
            fig = plt.figure()
            grph = fig.add_subplot(1,1,1)
            grph.set_autoscale_on(True)
            grph.set_title(r"$\chi^2$"+" vs. "+str(names[i]))
            grph.set_xlabel(str(names[i]))
            grph.set_ylabel(r"$\chi^2$")
            
            # finally, add the data
            grph.plot(vals[i],chisq,'ko')
            
            oname = str(conf.output)+"_"+str(names[i])+"_chisq.png"
            fig.savefig(oname,format="png")
            plt.close(fig)
        
    
    
def outputFitCurve(conf, d, m, mnum, chisq):
    fig, ax = plt.subplots()
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.set_autoscale_on(True)
    
    # set the plot title
    val = "%6.4f" %(chisq)
    ax.set_title('Fit of profile '+str(mnum)+' to '+str(conf.datafile)+r' $\chi^2$='+val)
    
    # set the text of the axes
    ax.set_xlabel('q  (1/${\AA}$)')
    ax.set_ylabel('Intensity  (1/cm)')
    
    # plot the data and the fit curve.
    ax.errorbar(d.x, d.y, yerr=d.dy, marker='o')
    # local_dy = np.zeros(len(m.x))
    # ax.errorbar(m.x, m.y, yerr=local_dy, color='r')
    ax.plot(m.x,m.y,color='r')
    
    oname = str(conf.output)+"%02d.png" %(mnum)
    fig.savefig(oname,format='png')
    plt.close(fig)
    
