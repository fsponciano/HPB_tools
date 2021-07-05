""" --- Energy level + absorption spectrum ('Big') diagram ---

	Written by J. Keaveney (2017)
	** Modified by F. S. Ponciano Ojeda (2021) **
	Quantum Light & Matter Group, Durham University, UK

	'Big' digram for showing the atomic energy levels, dipole-allowed transitions,
	and calculated absorption spectra at a given magnetic field strength.

	Calculation of atomic absorption spectra uses the *ElecSus* software, which
	can be found at https://github.com/jameskeaveney/ElecSus.
	Further details can be found in the following references:
	- Zentile, M. et al. ElecSus: A program to calculate the electric susceptibility
		of an atomic ensemble. Comp. Phys. Comm. 189 (2015), 162-174,
		http://dx.doi.org/10.1016/j.cpc.2014.11.023
	- Keaveney, J. et al. ElecSus: Extension to arbitrary geometry magneto-optics.
		Comp. Phys. Comm. 224 (2018), 311-324,
		https://doi.org/10.1016/j.cpc.2017.12.001
"""
import matplotlib

import numpy as np
import matplotlib.pyplot as plt
import time

# Fancy arrows
from matplotlib.patches import ConnectionPatch
import matplotlib.image as mpimg

# Parallel processing
from multiprocessing import Pool

import sys,os
from elecsus.libs.spectra import get_spectra
from elecsus.libs.spectra_energies import calc_chi_energies
import elecsus.libs.EigenSystem as ES

from scipy.constants import physical_constants, epsilon_0, hbar, c, e, h
kB = physical_constants['Boltzmann constant'][0]
mu_B = physical_constants['Bohr magneton in Hz/T'][0] / 1e9

# State decomposition code
from uncoupledbasis import AM_StateDecomp

from durhamcolours import *

# update matplotlib fonts etc
plt.rc('font',**{'family':'Serif','serif':['Times New Roman'],'weight':'bold'})
params={'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14,'legend.fontsize': 12,'mathtext.fontset':'cm','mathtext.rm':'serif'}
plt.rcParams.update(params)


def eval_energies(args):
	''' Use ElecSus to calculate the Hamiltonian for the atomic system and then
	diagonalise this to obtain the eigenenergies.'''

	isotope = args[0]
	Dline = args[1]
	Bfield = args[2]

	DlineHamiltonian = ES.Hamiltonian(isotope,Dline,1.0,Bfield)
	return DlineHamiltonian.groundEnergies, DlineHamiltonian.excitedEnergies

def big_diagram(BFIELD=1000):
	"""
		Main code to plot 'big' diagram with the following components:
		- Theoretical absorption spectrum (top panel)
		- Breit Rabi diagram for 0 to specified B-field (left)
		- Energy levels for ground and excited states (bottom panel)
		- Arrows for each transition, underneath the corresponding part of the spectrum
	"""

	##
	## First part - calculate the absorption spectrum
	##

	# Define the detuning axis based on what the magnetic field strength is (in GHz)
	Dmax = max(6,5 + (BFIELD/1e4 * 3 * mu_B))
	det_range = np.linspace(-Dmax,Dmax,int(3e4))

	# Input parameters to calculate the spectrum
	Bfield = BFIELD #alias
	ELEM = 'Rb'
	DLINE = 'D2'
	RB85FRAC = 0.0	# Pure Rb87
	LCELL = 1e-3
	TEMP = 100 # C ~ 373K

	# Voigt, horizontal polarisation
	pol = [1,0,0]
	p_dict = {'T':TEMP,'lcell':LCELL,'Elem':ELEM,'rb85frac':RB85FRAC,'Dline':DLINE,
		'Bfield':BFIELD,'Btheta':90*np.pi/180,'Bphi':45*np.pi/180,'BoltzmannFactor':True}

	[S0,S1,S2,S3] = get_spectra(det_range*1e3,pol,p_dict,outputs=['S0','S1','S2','S3'])

	lenergy87, lstrength87, ltransno87, lgl87, lel87, \
		renergy87, rstrength87, rtransno87, rgl87, rel87, \
		zenergy87, zstrength87, ztransno87, zgl87, zel87 = calc_chi_energies([1], p_dict)


	##
	## Second part - calculate the Breit-Rabi diagram
	##

	BreitRabiVals = np.linspace(0,BFIELD,2000)
	BreitRabiVals = np.append(BreitRabiVals,BreitRabiVals[-1])
	Bstep = BreitRabiVals[1] - BreitRabiVals[0]

	# Calculate Zeeman-shifted energy levels in parallel (uses multiprocessing module)
	po = Pool()
	res = po.map_async(eval_energies,(("Rb87","D2",BreitRabiVals[k],) for k in range(len(BreitRabiVals))))
	energies = res.get()
	gnd_energies = np.zeros((len(energies[0][0]),len(BreitRabiVals)))
	exc_energies = np.zeros((len(energies[0][1]),len(BreitRabiVals)))
	for jj, energyB in enumerate(energies):
		gnd_energies[:,jj] = energyB[0]
		exc_energies[:,jj] = energyB[1]
	po.close()
	po.join()

	# Energies at largest B-field value
	final_gnd_energies, final_exc_energies = eval_energies(("Rb87","D2",BreitRabiVals[-1]))


	##
	## Third part - calculate state decomposition
	##

	## Below values are for Rb-87. **Change for other atoms**.
	I=3.0/2; L=0; S=1.0/2; J=1.0/2
	output_states = AM_StateDecomp(I,L,S,J,atom='Rb',B=BFIELD/1e4)
	print('\nState decomposition at B = ',BFIELD/1e4)
	print(output_states)


	##
	## Fourth part - arrange the plot panels
	##

	fig = plt.figure("Big diagram at "+str(BFIELD/1e4)+' T',facecolor=None,figsize=(12,8))
	plt.clf()

	# Subplot arrangement
	xBR = 2
	xspec = 6
	yBRe = 3
	yBRg = 5
	yspec = 4

	xx = xBR + xspec
	yy = yBRe + yBRg + yspec

	ax_spec = plt.subplot2grid((yy,xx),(0,xBR),colspan=xspec,rowspan=yspec)
	ax_excBR = plt.subplot2grid((yy,xx),(yspec,0),colspan=xBR,rowspan=yBRe)
	ax_gndBR = plt.subplot2grid((yy,xx),(yspec+yBRe,0),colspan=xBR,rowspan=yBRg,sharex=ax_excBR)
	ax_eLev = plt.subplot2grid((yy,xx),(yspec,xBR),colspan=xspec,rowspan=yBRe,sharex=ax_spec,sharey=ax_excBR)
	ax_gLev = plt.subplot2grid((yy,xx),(yspec+yBRe,xBR),colspan=xspec,rowspan=yBRg,sharex=ax_spec,sharey=ax_gndBR)

	# Turn off axes for eLev and gLev axes
	for ax in [ax_eLev,ax_gLev]:
		ax.set_frame_on(False)
		for parameter in [ax.get_xticklabels(),ax.get_yticklabels(),ax.get_xticklines(),ax.get_yticklines()]:
			plt.setp(parameter,visible=False)

	plt.setp(ax_excBR.get_xticklabels(),visible=False)

	ax_excBR.spines['right'].set_color('none')
	ax_gndBR.spines['right'].set_color('none')
	ax_gndBR.spines['top'].set_color('none')
	ax_excBR.spines['top'].set_color('none')
	ax_excBR.spines['bottom'].set_color('none')
	ax_gndBR.xaxis.set_ticks_position('bottom')
	ax_excBR.xaxis.set_ticks_position('none')

	ax_excBR.tick_params(axis='y',left=True,right=False)
	ax_gndBR.tick_params(axis='y',left=True,right=False)

	# axis labels
	ax_spec.set_xlabel('Detuning (GHz)')
	ax_spec.xaxis.set_label_position('top')
	ax_spec.tick_params(axis='x',bottom=True,top=True,labelbottom=False,labeltop=True)

	ax_spec.set_ylabel('Transmission')
	ax_excBR.set_ylabel('$5P_{3/2}$ energy (GHz)')
	ax_gndBR.set_ylabel('$5S_{1/2}$ energy (GHz)')

	ax_gndBR.set_xlabel('Magnetic Field (T)')

	fig.subplots_adjust(left=0.07,right=0.98,top=0.93,bottom=0.085,hspace=0.34,wspace=0)

	#Ghost axes for actually plotting the Breit-Rabi data
	eleft = ax_excBR.get_position().extents[0:2]
	eright = ax_eLev.get_position().extents[2:]
	gleft = ax_gndBR.get_position().extents[0:2]
	gright = ax_gLev.get_position().extents[2:]

	ax_e_bound = np.append(eleft,eright-eleft)
	ax_g_bound = np.append(gleft,gright-gleft)

	print('\nAxes bounds for B-R diagram:')
	print(ax_e_bound)
	print(ax_g_bound)

	ax_e = fig.add_axes(ax_e_bound,frameon=False,facecolor=None)
	ax_g = fig.add_axes(ax_g_bound,frameon=False,facecolor=None)

	ax_g.set_xticks([])
	ax_g.set_yticks([])
	ax_e.set_xticks([])
	ax_e.set_yticks([])


	##
	## Fifth part - Add the data to the figure
	##

	# Edit last magnetic field value
	BreitRabiVals[-1] = BreitRabiVals[-2] * ((xspec + xBR) / xBR)
	print('\nMagnetic field values (Breit-Rabi diagram)')
	print(BreitRabiVals)
	ax_spec.plot(det_range,S0.real,lw=2,color=d_black)

	#convert to GHz from MHz
	exc_energies /= 1e3
	gnd_energies /= 1e3
	final_exc_energies /= 1e3
	final_gnd_energies /= 1e3

	for energy in exc_energies[int(len(final_exc_energies)/3):]:
		ax_e.plot(BreitRabiVals/1e4,energy,color=d_black,lw=1)
	for energy in gnd_energies:
		ax_g.plot(BreitRabiVals/1e4,energy,color=d_black,lw=1.5)

	ax_excBR.set_xlim(0,(Bfield + 10*Bstep)/1e4)
	for ax in [ax_g,ax_e]:
		ax.set_ylim(ax.get_ylim()[0]*1.15,ax.get_ylim()[1]*1.15)
		ax.set_xlim(BreitRabiVals[0]/1e4, BreitRabiVals[-1]/1e4)
	ax_excBR.set_ylim(ax_e.get_ylim())
	ax_gndBR.set_ylim(ax_g.get_ylim())

	ax_spec.set_xlim(det_range[0],det_range[-1])
	ax_spec.set_ylim(ax_spec.get_ylim()[0],1.01)


	##
	## Sixth part - Add arrows for each transition
	##

	print('Sigma minus transitions:')
	print(sorted(lenergy87))
	print('Sigma plus transitions:')
	print(sorted(renergy87))
	print('Pi transitions:')
	print(sorted(zenergy87))

	for energy in lenergy87:
		ax_spec.axvline(energy/1e3,color=d_purple,lw=1.5)
	for energy in renergy87:
		ax_spec.axvline(energy/1e3,color=d_blue,lw=1.5)
	for energy in zenergy87:
		ax_spec.axvline(energy/1e3,color=d_olive,lw=1.5,linestyle='dashed')

	# Coordinates for arrows - sigma minus transitions (purple)
	xy1s = zip(lenergy87/1e3,lgl87/1e3)
	xy2s = zip(lenergy87/1e3,lel87/1e3)
	ecol = d_purple
	fcol = 0.5 * (np.array(d_lightpurple) + np.array(d_purple))
	alpha = 0.9
	#styles = ['solid','solid','solid','solid','dashed','dashed','dashed','dashed']
	for xy1,xy2,strength in zip(xy1s,xy2s,lstrength87):
		#if (xy1[0] > 15) or (xy1[0]<-15):
		coordsA = 'data'
		coordsB = 'data'
		con = ConnectionPatch(xy1,xy2,coordsA,coordsB,
					arrowstyle="simple",shrinkB=0,
					axesA=ax_gLev,axesB=ax_eLev,mutation_scale=25,
					ec=ecol,fc=fcol,lw=1.25,alpha=alpha)
		ax_gLev.add_artist(con)

	# Coordinates for arrows - sigma plus transitions (blue)
	xy1s = zip(renergy87/1e3,rgl87/1e3)
	xy2s = zip(renergy87/1e3,rel87/1e3)
	ecol = d_blue
	fcol = 0.5 * (np.array(d_midblue) + np.array(d_blue))
	alpha = 0.9
	#styles = ['solid','solid','solid','solid','dashed','dashed','dashed','dashed']
	for xy1,xy2,strength in zip(xy1s,xy2s,rstrength87):
		#if (xy1[0] > 15) or (xy1[0]<-15):
		coordsA = 'data'
		coordsB = 'data'
		con = ConnectionPatch(xy1,xy2,coordsA,coordsB,
					arrowstyle="simple",shrinkB=0,
					axesA=ax_gLev,axesB=ax_eLev,mutation_scale=25,
					ec=ecol,fc=fcol,lw=1.25,alpha=alpha)
		ax_gLev.add_artist(con)

	# Coordinates for arrows - pi transitions (olive)
	xy1s = zip(zenergy87/1e3,zgl87/1e3)
	xy2s = zip(zenergy87/1e3,zel87/1e3)
	ecol = d_darkolive
	fcol = d_olive#darkyellow#olive #(0.16,0.85,0.16)
	alpha = 0.6
	#styles = ['solid','solid','solid','solid','dashed','dashed','dashed','dashed']
	for xy1,xy2,strength in zip(xy1s,xy2s,zstrength87):
		#if (xy1[0] < 15) and (xy1[0]>-15):
		coordsA = 'data'
		coordsB = 'data'
		con = ConnectionPatch(xy1,xy2,coordsA,coordsB,
					arrowstyle="simple",shrinkB=0,
					axesA=ax_gLev,axesB=ax_eLev,mutation_scale=25,
					ec=ecol,fc=fcol,lw=1.25,alpha=alpha)
		ax_gLev.add_artist(con)

	# Add B-field info to plot - top left
	fig.text(0.1,0.78-0.03,'L = '+str(LCELL*1e3)+' mm',size=18,ha='center')
	fig.text(0.1,0.82-0.03,r'T = '+str(TEMP)+' $^{\circ}$C',size=18,ha='center')
	fig.text(0.1,0.86-0.03,'B = '+str(Bfield/1e4)+' T',size=18,ha='center')
	fig.text(0.1,0.90-0.03,str(DLINE)+' Line',size=18,ha='center')
	fig.text(0.1,0.94-0.03,'$^{87}$Rb',size=18,ha='center')


	##
	## Finally - show the plot and save the figure
	##

	ax_spec.set_xlim(-Dmax,Dmax)

	# fig.savefig('./BR_plotoutput'+str(Bfield)+'.pdf',dpi=300)
	# fig.savefig('./BR_plotoutput'+str(Bfield)+'.png',dpi=300)

	plt.show()

	print('--- End of calculations ---')

def make_vid_stills():
	Brange = np.arange(0,3000,200)
	Brange[0] = 1
	Brange = np.append(Brange,np.arange(3000,15000,1000))
	Brange = np.append(Brange,np.arange(15000,85000,5000))
	print('Number of frames: ',len(Brange))

	for B in Brange:
		big_diagram(B)

def state_decomp(Bfield,T=120,L=1,Pol=50):
	""" Make Briet-Rabi diagram with Boltzmann factor """

	# Calculate Boltzmann factor
	p = ['Rb','D2',Bfield,T,1,0,T,45,Pol,0,30,True,0,0,True]
	bd, ge = calculate(np.array([1]), p, OutputType='Boltz')

	# Calculate B-R diagram
	BreitRabiVals = np.linspace(0,Bfield,2000)
	BreitRabiVals = np.append(BreitRabiVals,BreitRabiVals[-1])
	Bstep = BreitRabiVals[1] - BreitRabiVals[0]

	# do it in parallel! (Multiprocessing module)
	po = Pool()
	res = po.map_async(eval_energies,(("Rb87","D2",BreitRabiVals[k],) for k in range(len(BreitRabiVals))))
	energies = res.get()
	g_energies = np.zeros((len(energies[0][0]),len(BreitRabiVals)))
	e_energies = np.zeros((len(energies[0][1]),len(BreitRabiVals)))
	for jj, energyB in enumerate(energies):
		g_energies[:,jj] = energyB[0]
		e_energies[:,jj] = energyB[1]
	po.close()
	po.join()

	## calculate state decomposition

	## Rb-87 only!
	I = 3./2; L = 0; S=1./2; J=1./2
	output_states = AM_StateDecomp(I,L,S,J,atom='Rb',B=Bfield/1e4)
	print(output_states)


	## make figure
	fig = plt.figure("State decomposition factor",figsize=(7,7))
	fig.clf()

	ax = fig.add_subplot(111)

	fig.subplots_adjust(left=0.12,right=0.95,top=0.93,bottom=0.08,hspace=0.34,wspace=0)

	ax.set_xlabel('Magnetic field (G)')
	ax.set_ylabel(r'$5S_{1/2}$ energy (GHz)')

	for en in g_energies:
		ax.plot(BreitRabiVals,en/1e3,'k',lw=1.5)

	for i in range(len(output_states)):
		energy = output_states[i][0]
		# look for states with same energy

		if  i != len(output_states)-1:
			if output_states[i][0] == output_states[i+1][0]:
				# concatenate strings
				textstr = str(round(float(output_states[i][1]),2))+' '+output_states[i][2]+' + ' \
					+str(round(float(output_states[i+1][1]),2))+' '+output_states[i+1][2]
			else:
				textstr = output_states[i][2]

		if output_states[i][0] == output_states[i-1][0]:
			textstr = ''

		ax.text(Bfield*1.05,energy,textstr,size=13,va='center')

	ax.set_xlim(0,Bfield*2)

	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')

	ax.tick_params(axis='x',bottom=True,top=False)
	ax.tick_params(axis='y',left=True,right=False)

	#fig.savefig('./movie_stills_Bfactor/BoltzmannPlot'+str(Bfield)+'.pdf')
	fig.savefig('./movie_stills_decomposition/BoltzmannPlot'+str(Bfield)+'.png')

def make_decomp_vid_stills():
	Brange = np.arange(0,3000,200)
	Brange[0] = 10
	Brange = np.append(Brange,np.arange(3000,15000,1000))
	Brange = np.append(Brange,np.arange(15000,85000,5000))
	print('Number of frames:', len(Brange), '\n\n\n')

	for B in Brange:
		state_decomp(B)

def Boltz_factor(Bfield,T=120, L=1, Pol=50):
	""" Make Breit-Rabi diagram with Boltzmann factor """

	## Calculate BR

	BreitRabiVals = np.linspace(0,Bfield,2000)
	BreitRabiVals = np.append(BreitRabiVals,BreitRabiVals[-1])
	Bstep = BreitRabiVals[1] - BreitRabiVals[0]

	# do it in parallel! (Multiprocessing module)
	po = Pool()
	res = po.map_async(eval_energies, (("Rb87","D2",BreitRabiVals[k],) for k in range(len(BreitRabiVals))))
	energies = res.get()
	g_energies = np.zeros((len(energies[0][0]),len(BreitRabiVals)))
	e_energies = np.zeros((len(energies[0][1]),len(BreitRabiVals)))
	for jj, energyB in enumerate(energies):
		g_energies[:,jj] = energyB[0]
		e_energies[:,jj] = energyB[1]
	po.close()
	po.join()


	## calculate state decomposition
	## Rb-87 only!
	I = 3./2; L = 0; S=1./2; J=1./2
	output_states = AM_StateDecomp(I,L,S,J,atom='Rb',B=Bfield/1e4)
	print(output_states)

	#calculate Boltzmann factor
	p = ['Rb','D2',Bfield,T,1,0,T,45,Pol,0,30,True,0,0,True]
	bd, ge = calculate(np.array([1]), p, OutputType='Boltz')

	print(bd, ge)
	bd = np.array(bd)
	bd /= bd.max()

	ge = np.array(ge)/1e3 + g_energies.min()/1e3

	## make figure
	fig = plt.figure("Boltzmann factor",figsize=(7,7))
	fig.clf()

	ax = fig.add_subplot(111)

	fig.subplots_adjust(left=0.12,right=0.95,top=0.93,bottom=0.08,hspace=0.34,wspace=0)

	ax.set_xlabel('Magnetic field (G)')
	ax.set_ylabel(r'$5S_{1/2}$ energy (GHz)')

	for en in g_energies:
		ax.plot(BreitRabiVals,en/1e3,'k',lw=1.5)

	for b,e in zip(bd,ge):
		ax.text(Bfield*1.025,e,str(round(b,4)),size=13,va='center')

	ax.text(Bfield*1.025,ax.get_ylim()[1],'Boltzmann Factor')

	fig.text(0.15,0.9,'B = '+str(Bfield/1e4)+' T')
	ax.set_xlim(0,Bfield*1.35)

	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')

	ax.tick_params(axis='x',bottom=True,top=False)
	ax.tick_params(axis='y',left=True,right=False)

	#fig.savefig('./movie_stills_Bfactor/BoltzmannPlot'+str(Bfield)+'.pdf')
	fig.savefig('./movie_stills_Bfactor/BoltzmannPlot'+str(Bfield)+'.png')

def make_Bfactor_stills():
	Brange = np.arange(0,3000,200)
	Brange[0] = 10
	Brange = np.append(Brange,np.arange(3000,15000,1000))
	Brange = np.append(Brange,np.arange(15000,85000,5000))
	print('Number of frames:', len(Brange), '\n\n\n')

	for B in Brange:
		Boltz_factor(B)

if __name__ == '__main__':
	Bfields = [4000,15000]
	for Bfield in Bfields:
		big_diagram(Bfield)
	#make_vid_stills()
