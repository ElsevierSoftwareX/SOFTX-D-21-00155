r"""
sas_temper

sa_calc.py:  this code is for the intensity and	chi-squared calculations.  
             The smearing is also handled here.  This code will work with 
             sasmodels and appropraitely-defined custom models.  There is 
             no GPU parallelization here.  The smearing is done via a simple
             integral, but the code will be expanded to use the NIST
             Gaussian integration routine built into sasmodels.
             
Oak Ridge National Laboratory, 2020

"""

import numpy as np
import math as m
from sasmodels.sasview_model import _make_standard_model
from sasmodels.weights import *

import sas_temper.sas_data as sas_data
import sas_temper.modelconfig as modelconfig

# d is the input data, m is the model configuration and usm is the unsmeared profile
def calc_profile(d, mc, usm):
	# this function is only used to calculate 
	# a smeared profile from an unsmeared one
	local = sas_data.Model(d, unsmeared = False)
	
	for i in range(0,len(d.x)):
		tot = 0.0
		min = d.x[i] - 4.5*d.dx[i]
		qstep = 4.5*d.dx[i]/30.0
		for j in range(0,61):
			qv = m.fabs(min + j*qstep)
			tmp = np.interp(qv,usm.x,usm.y)
			# the number on the next line is ~sqrt(2*Pi)
			tot += tmp*qstep*m.exp(-(qv-d.x[i])*(qv-d.x[i])/(2.0*d.dx[i]*d.dx[i]))/(2.506628273*d.dx[i])
			
		local.y[i] = tot
		local.x[i] = d.x[i]
		
	return local

# d is the input data, m is the model configuration	
def calc_profile_usm(d, mc):
	# configure and compute the model form factor
	local = sas_data.Model(d, unsmeared = False)
	
	# this local configuration is where we set the log and coupled parameters
	# into values that can be used for the intensity calculation
	lc = modelconfig.convert_conf(mc)
	
	Profile = _make_standard_model(mc.name)
	model = Profile()
	for i, p in enumerate(lc.params):
		model.setParam(p.name,p.val)
		
		if p.polydispersity is not None:
			if p.polydispersity=="SchulzDispersion":
				model.set_dispersion(p.name,SchulzDispersion(npts=70,width=p.polydispersity.val,nsigmas=15))
			elif p.polydispersity=="GaussianDispersion":
				model.set_dispersion(p.name,GaussianDispersion(npts=35,width=p.polydispersity.val,nsigmas=3))
			elif p.polydispersity=="LogNormalDispersion":
				model.set_dispersion(p.name,LogNormalDispersion(npts=35,width=p.polydispersity.val,nsigmas=3))
			elif p.polydispersity=="RectangleDispersion":
				model.set_dispersion(p.name,RectangleDispersion(npts=35,width=p.polydispersity.val,nsigmas=2))
			else:
				#default to Gaussian
				model.set_dispersion(p.name,GaussianDispersion(npts=35,width=p.polydispersity.val,nsigmas=3))
		
	# this line actually does the form factor calculation
	local.y = model.evalDistribution(local.x)
	
	# configure and compute the structure factor if we have one
	if lc.sq is not None:
		local_sq = sas_data.Model(d, unsmeared = False)
		
		SQProfile = _make_standard_model(lc.sq.type)
		sqmodel = SQProfile()
		for i,sqp in enumerate(lc.sq.params):
			sqmodel.setParam(sqp.name,sqp.val)
			
			if sqp.polydispersity is not None:
				if sqp.polydispersity=="SchulzDispersion":
					model.set_dispersion(sqp.name,SchulzDispersion(npts=70,width=sqp.polydispersity.val,nsigmas=15))
				elif sqp.polydispersity=="GaussianDispersion":
					model.set_dispersion(sqp.name,GaussianDispersion(npts=35,width=sqp.polydispersity.val,nsigmas=3))
				elif sqp.polydispersity=="LogNormalDispersion":
					model.set_dispersion(sqp.name,LogNormalDispersion(npts=35,width=sqp.polydispersity.val,nsigmas=3))
				elif sqp.polydispersity=="RectangleDispersion":
					model.set_dispersion(sqp.name,RectangleDispersion(npts=35,width=sqp.polydispersity.val,nsigmas=2))
				else:
					#default to Gaussian
					model.set_dispersion(sqp.name,GaussianDispersion(npts=35,width=sqp.polydispersity.val,nsigmas=3))
					
		# this line actually does the structure factor calculation
		local_sq.y = sqmodel.evalDistribution(local_sq.x)
		
		# and this makes P(q)*S(q)
		local.y = local.y*local_sq.y
		
	# and out go the results
	return local
	
# mc is the model configuration, d is the input data, mi is the calculated model profile
def chisq(mc, d, mi):
	# this is quite straightforward
	sum = 0
	for i in range(0,len(d.x)):
		sum += (d.y[i]-mi.y[i])*(d.y[i]-mi.y[i])/(d.dy[i]*d.dy[i])
		
	sum = m.sqrt(sum/(mc.free_params-1))
	
	return sum
	
		
