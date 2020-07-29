r"""
sas_temper

sa_data.py:  this is a container class for the data that allows 
             the program to drop data outside the q-range of 
             interest, which the sasmodels does not directly support
             in a clear manner.

Oak Ridge National Laboratory, 2020

"""

import numpy as np
import math as m

class SAData(object):
	def __init__(self, filename, qmin=0.0, qmax=1.0):
		
		self.name = str(filename)
		
		# set up our q-range for the modeling
		self.qmin = float(qmin)
		self.qmax = float(qmax)
		
		#print("Input file: "+self.name+"; qmin "+str(self.qmin)+"; qmax "+str(self.qmax))
		
		with open(self.name, 'r') as f:
			flines = f.readlines()
			
			data_start = 0
			for line in flines:
				if line.startswith("#"):
					data_start += 1
				else:
					break
					
			split_line = flines[data_start].split()
			columns = len(split_line)
			rows = len(flines) - data_start
			#print("The number of columns is "+str(columns)+" and the number of rows is "+str(rows))
			#print("data_start "+str(data_start)+"; lines in the file "+str(len(flines)))
			
			dvals = 0
			for i in range(data_start,len(flines)):
				line = flines[i].strip().strip('\n')
				
				if line:
					#sline = flines[i].split()
					sline = line.split()
					
					if float(sline[0]) >= self.qmin :
						if float(sline[0]) <= self.qmax :
							dvals += 1
                            
			self.x = np.zeros(dvals)
			self.y = np.zeros(dvals)
			if columns is 2:
				self.dy = None
				self.dx = None
			elif columns is 3:
				self.dy = np.zeros(dvals)
				self.dx = None
			elif columns is 4:
				self.dy = np.zeros(dvals)
				self.dx = np.zeros(dvals)
			
			j = 0
			for i in range(data_start,len(flines)):
				line = flines[i].strip().strip('\n')
				
				if line:
					#sline = flines[i].split()
					sline = line.split()
					
					if float(sline[0]) >= self.qmin :
						if float(sline[0]) <= self.qmax :
							self.x[j] = float(sline[0])
							self.y[j] = float(sline[1])
							if columns is 3:
								self.dy[j] = float(sline[2])
							elif columns is 4:
								self.dy[j] = float(sline[2])
								self.dx[j] = float(sline[3])
							j += 1
				
		f.close()

# This is a model intensity profile that can be smeared or unsmeared
class Model(object):
	# construct the object depending on whether 
	# it contains a smeared or unsmeared profile
	def __init__(self, data, unsmeared=False):
		if unsmeared:
			qvals = 10*len(data.x)
			qmin = -4.0
			qmax = m.log10(2.0*np.amax(data.x))
			self.x = np.logspace(qmin,qmax,qvals,base=10.0,dtype=np.float64)
			self.y = np.zeros(qvals)
		else:
			self.x = np.zeros(len(data.x))
			self.x = data.x
			self.y = np.zeros(len(data.x))
			
