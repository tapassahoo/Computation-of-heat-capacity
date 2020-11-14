import os                                                                                                                       
import matplotlib
import matplotlib.axes
# matplotlib.use('eps')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.ticker as mticker
import numpy as np
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
from pylab import *

rc('text', usetex=True)
size=24
params = {'legend.fontsize': size*0.6,
    'figure.figsize': (8,6),
    'axes.labelsize': size,
    'axes.titlesize': size,
    'xtick.labelsize': size*0.75,
    'ytick.labelsize': size*0.75,
    'axes.titlepad': size}
plt.rcParams.update(params)
matplotlib.rcParams.update({'font.size': size*0.75})                                                                            
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

valmin=0.1
valmax=100
dt = 0.1
nt = int(((valmax-valmin)+dt*0.5)/dt)
nt+=1
tempList = [valmin+dt*i for i in range(nt)]

colorList = ['r', 'm', 'g', 'b']
markerList = ['o', 's', 'p', '<', '<', '8', 'p']
lsList = ['-', '--', '-.', ':']

rptList = [3.0]
interaction = 'cost-1'
#interaction = 'free'
output_prefix = '/Users/tsahoo/ResultsOfExact/'
if (interaction == 'free'):
	File_prefix = '/Users/tsahoo/ResultsOfExact/eigen-values-of-1-'
	output_prefix += 'Cv-of-free-'
else:
	File_prefix = '/Users/tsahoo/ResultsOfExact/eigen-values-of-2-'
	output_prefix += 'Cv-of-'+interaction+'-'

isomerList = ["equilibrium", "spinless", "para", "ortho"]
for pltid, isomer in enumerate(isomerList):

	if (isomer == "spinless"):
		spin=""
		spin_label="Spinless-"
	if (isomer == "para"):
		spin="p-"
		spin_label="Para-"
	if (isomer == "ortho"):
		spin="o-"
		spin_label="Ortho-"
	if (isomer == "equilibrium"):
		spin_label="Equilibrium-"

	for rcom1 in rptList:
		rcom="{:3.2f}".format(rcom1)

		if (interaction == 'free'):
			File_para = File_prefix+'p-H2O-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
		else:
			File_para = File_prefix+'p-H2O-one-rotor-fixed-'+interaction+'-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
		valTotalEnergy_para = np.genfromtxt(File_para, unpack=True, usecols=[0], skip_header=0, skip_footer=0)

		if (isomer != "equilibrium"):
			if (interaction == 'free'):
				FileToBeRead = File_prefix+spin+'H2O-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
			else:
				FileToBeRead = File_prefix+spin+'H2O-one-rotor-fixed-'+interaction+'-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
			valTotalEnergy = np.genfromtxt(FileToBeRead, unpack=True, usecols=[0], skip_header=0, skip_footer=0)
			print(valTotalEnergy_para[0])
			valTotalEnergy = valTotalEnergy-valTotalEnergy_para[0]
			valTotalEnergySq=np.square(valTotalEnergy)

			cvList = np.zeros(nt,float)
			for i, temperature in enumerate(tempList):
				beta=1.0/temperature
				exp1=np.exp(-beta*valTotalEnergy)
				
				avgEnergy=np.dot(valTotalEnergy,exp1)/np.sum(exp1)
				avgEnergySq=np.dot(valTotalEnergySq,exp1)/np.sum(exp1)
				cvList[i]=beta*beta*(avgEnergySq-avgEnergy*avgEnergy)
				tempList[i]=temperature
				
			cv_write=np.array([tempList,cvList])
			filename = output_prefix+spin+'H2O-Rpt'+rcom+'Angstrom.txt'
			np.savetxt(filename,cv_write.T,fmt='%20.8f',delimiter=' ',header='First col: Temperature in Kelvin; second col: Cv')
			
			labelstr = spin_label+r'$\mathrm{H_2 O}$'
			plt.plot(tempList, cvList, color=colorList[pltid], ls=lsList[pltid], linewidth=1,  marker=markerList[pltid], markersize=4, label=labelstr)

		if (isomer == "equilibrium"):
			if (interaction == 'free'):
				File_para = File_prefix+'p-H2O-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
			else:
				File_para = File_prefix+'p-H2O-one-rotor-fixed-'+interaction+'-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
			valTotalEnergy_para = np.genfromtxt(File_para, unpack=True, usecols=[0], skip_header=0, skip_footer=0)
			val_para = valTotalEnergy_para-valTotalEnergy_para[0]
			valTotalEnergySq_para = np.square(val_para)

			if (interaction == 'free'):
				File_ortho = File_prefix+'o-H2O-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
			else:
				File_ortho = File_prefix+'o-H2O-one-rotor-fixed-'+interaction+'-jmax12-Rpt'+rcom+'Angstrom-grids-27-50-diag.txt'
			valTotalEnergy_ortho = np.genfromtxt(File_ortho, unpack=True, usecols=[0], skip_header=0, skip_footer=0)
			val_ortho = valTotalEnergy_ortho-valTotalEnergy_para[0]
			valTotalEnergySq_ortho = np.square(val_ortho)

			cvList = np.zeros(nt)
			for i, temperature in enumerate(tempList):
				beta=1.0/temperature
				exp1_ortho=np.exp(-beta*val_ortho)
				exp1_para=np.exp(-beta*val_para)
				partition_fun = 3*np.sum(exp1_ortho)+np.sum(exp1_para)
				
				avgEnergy=(3*np.dot(val_ortho,exp1_ortho)+np.dot(val_para,exp1_para))/partition_fun
				avgEnergySq=(3*np.dot(valTotalEnergySq_ortho,exp1_ortho)+np.dot(valTotalEnergySq_para,exp1_para))/partition_fun
				cvList[i]=beta*beta*(avgEnergySq-avgEnergy*avgEnergy)
				tempList[i]=temperature
				
			cv_write=np.array([tempList,cvList])
			filename = output_prefix+'equilibrium-H2O-Rpt'+rcom+'Angstrom.txt'
			np.savetxt(filename,cv_write.T,fmt='%20.8f',delimiter=' ',header='First col: Temperature in Kelvin; second col: Cv')
			
			labelstr = spin_label+r'$\mathrm{H_2 O}$'
			plt.plot(tempList, cvList, color=colorList[pltid], ls=lsList[pltid], linewidth=1,  marker=markerList[pltid], markersize=4, label=labelstr)



plt.xlabel(r'$\mathrm{Temperature} \ (\mathrm{Kelvin})$', labelpad=5)            
plt.ylabel(r'$\mathrm{C_{V}} \ (\mathrm{k_{B}})$',labelpad=5)


labelListPIGS = {0:'PIGS: ('+r'$\pi$'+' 0, 0)', 1:'PIGS: ('+'0, 0, 0)'}

#For text lebelling
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
midpointy = 0.5*(ymax-ymin)
deltay = midpointy*0.15
xmin, xmax = plt.xlim()
midpointx = 0.5*(xmax-xmin)
deltax = midpointx*0.15
textpositionx = xmin+midpointx-0.25*midpointx
textpositiony = ymin+midpointy
																															
rcom1="{:3.1f}".format(rcom1)
if (interaction == 'free'):
	plt.text(xmin+(xmax-xmin)*0.3,ymax-(ymax-ymin)*0.1,r'$N = \ 1$; No interaction')
elif (interaction == 'cost-1'):
	#plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.5,r'$N = \ 1$; '+r'$(\pi$'+' 0, 0)')
	#plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.6,r'$r = $'+rcom1+r'$\mathrm{\AA}$')
	plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.1,r'$N = \ 1$; '+r'$(\pi$'+' 0, 0)')
	plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.2,r'$r = $'+rcom1+r'$\mathrm{\AA}$')
elif (interaction == 'cost1'):
	#plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.5,r'$N = \ 1$; (0, 0, 0)')
	#plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.6,r'$r = $'+rcom1+r'$\mathrm{\AA}$')
	plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.1,r'$N = \ 1$; (0, 0, 0)')
	plt.text(xmin+(xmax-xmin)*0.4,ymax-(ymax-ymin)*0.2,r'$r = $'+rcom1+r'$\mathrm{\AA}$')

#For ticks manipulating                                                                                                     
plt.minorticks_on()
plt.tick_params(axis="both", direction="in", which="minor", right=False, top=False, length=2)
plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)
#Adjust the plot
plt.legend(numpoints=1,loc=('center right'))
#plt.legend(numpoints=1,loc=('upper right'))
plt.subplots_adjust(top=0.99,bottom=0.12,left=0.11,right=0.99,hspace=0.0,wspace=0.0)

FilePlot=output_prefix+'H2O-Rpt'+rcom+'Angstrom'
plt.savefig(FilePlot+".pdf", format='pdf') 
plt.show()
