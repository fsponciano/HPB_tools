"""

BIG DIAGRAM for All transitions (including pi transitions) at 1.5 T
For Boltzmann proposal paper

"""

import matplotlib
#matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import time

#animation (needed?)
import matplotlib.animation as ani

# parallel processing
from multiprocessing import Pool

import sys
# sys.path.append('E:\James\Documents\Programming\ElecSus arb Bfield\elecsus\libs')
sys.path.append('C:/Users/Francisco S Ponciano/Downloads/ElecSus-master_3.0.6_analytic/ElecSus-master/elecsus/')
from libs.spectra import get_spectra
from libs.spectra_energies import calc_chi_energies
#sys.path.append('E:\James\Documents\Programming\ElecSus Boltzmann\elecsus\libs')
import libs.EigenSystem as ES
#import EigenSystem_energies as ES

#fancy arrows
from matplotlib.patches import ConnectionPatch
import matplotlib.image as mpimg


# elecsus - Boltzmann version (local, not the installed version)
#from elecsus_methods import calculate
#from elecsus_methods_energies import calculate_energies
#from libs import EigenSystem_energies as ES
from scipy.constants import physical_constants, epsilon_0, hbar, c, e, h
kB = physical_constants['Boltzmann constant'][0]

# Dan Whiting's state decomposition code
from uncoupledbasis import MomentumDecomp

from durhamcolours import *

import time
import os

# update matplotlib fonts etc
plt.rc('font',**{'family':'Serif','serif':['Times New Roman']})
params={'axes.labelsize':13,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize': 11}
plt.rcParams.update(params)

def main():
	d = np.linspace(-80,80,10000)*1e3
	d_zoom = np.linspace(-80,80,10000)

	BFIELD = 20000
	TEMP = 100
	LCELL = 1e-3

	ELEM = 'Rb'
	DLINE = 'D2'
	RB85FRAC = 1.0 # Nat abundance?

	# Voigt, 45 degree linear pol

	#Boltzmann ON
	pol = 1./np.sqrt(2)*np.array([1.0,-1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM,'BoltzmannFactor':True}
	[S0_on,S3_on] = get_spectra(d,pol,p_dict,outputs=['S0','S3'])

	#Boltzmann OFF
	pol = 1./np.sqrt(2)*np.array([1.0,-1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM,'BoltzmannFactor':False}
	[S0_off,S3_off] = get_spectra(d,pol,p_dict,outputs=['S0','S3'])


	#Boltzmann ON
	pol = 1./np.sqrt(2)*np.array([1.0,-1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM,'BoltzmannFactor':True}
	[S0_on_zoom,S3_on_zoom] = get_spectra(d_zoom,pol,p_dict,outputs=['S0','S3'])

	#Boltzmann OFF
	pol = 1./np.sqrt(2)*np.array([1.0,-1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM,'BoltzmannFactor':False}
	[S0_off_zoom,S3_off_zoom] = get_spectra(d_zoom,pol,p_dict,outputs=['S0','S3'])

	## same again for zoom plots
	##
	##
	##


	'''
	fig1 = plt.figure(figsize=(7,9))
	ax1 = fig1.add_subplot(411)
	ax2 = fig1.add_subplot(412,sharex=ax1)
	ax3 = fig1.add_subplot(413,sharex=ax1)
	ax4 = fig1.add_subplot(414,sharex=ax1)

	ax1.plot(d/1e3,S0_on,color=d_purple,lw=2)
	ax2.plot(d/1e3,S1_on,color=d_purple,lw=2)
	ax3.plot(d/1e3,S2_on,color=d_purple,lw=2)
	ax4.plot(d/1e3,S3_on,color=d_purple,lw=2)

	ax1.set_xlim(-80,80)

	ax1.set_ylim(0.45,1.02)
	ax2.set_ylim(-0.55,0.55)
	ax3.set_ylim(-1.02,1.02)
	ax4.set_ylim(-1.02,1.02)

	ax4.set_xlabel('Detuning (GHz)')
	ax1.set_ylabel('S0')
	ax2.set_ylabel('S1')
	ax3.set_ylabel('S2')
	ax4.set_ylabel('S3')

	plt.tight_layout()
	'''

	fig2 = plt.figure(figsize=(5,7*5./6))

	N = 3
	gap = 0
	yy = 3*N + gap
	xx = 7

	ax1a = plt.subplot2grid((yy,xx),(0,0),rowspan=N,colspan = xx)
	ax2a = plt.subplot2grid((yy,xx),(N,0),rowspan=N,colspan = xx,sharex=ax1a)
	ax3a = plt.subplot2grid((yy,xx),(2*N+gap,1),rowspan=N,colspan = xx-2)
	#fig2.add_subplot(311)
	#ax2a = fig2.add_subplot(312,sharex=ax1a)
	#ax3a = fig2.add_subplot(313)

	fig2.subplots_adjust(bottom=0.085,top=0.93,right=0.96,left=0.11,hspace=0.36)


	ax1a.plot(d/1e3,S0_on,color=d_purple,lw=2,label=r'Thermally populated')
	ax1a.plot(d/1e3,S0_off,'--',color=d_olive,lw=2,label=r'Equally populated')

	ax2a.plot(d/1e3,S3_on,color=d_purple,lw=2,label=r'ON')
	ax2a.plot(d/1e3,S3_off,'--',color=d_olive,lw=2,label=r'OFF')

	ax3a.plot(d_zoom,S3_on_zoom,color=d_purple,lw=3,label=r'Thermally populated')
	ax3a.plot(d_zoom,S3_off_zoom,'--',color=d_olive,lw=3,label=r'Equally populated')

	ax2a.axvspan(-0.2,0.2,color=d_midblue,alpha=0.7)

	ax1a.xaxis.set_label_position('top')
	ax1a.tick_params(axis='x',bottom=True,top=True,labelbottom=False,labeltop=True)
	plt.setp(ax2a.get_xticklabels(),visible=False)

	ax1a.set_xlabel('Detuning (GHz)')
	ax3a.set_xlabel('Detuning (MHz)')

	ax1a.set_ylabel('Transmission, S0')
	ax2a.set_ylabel('S3')
	labelx=-0.075
	ax1a.yaxis.set_label_coords(labelx, 0.5)
	ax2a.yaxis.set_label_coords(labelx, 0.5)


	ax2a.set_xlim(d[0]/1e3,d[-1]/1e3)
	ax3a.legend(bbox_to_anchor=(1.25, 0.98))

	ax1a.set_xlim(-15,15)
	ax1a.set_ylim(0.5,1.02)

	ax2a.axhline(0,color='k',linestyle='dashed')
	ax3a.axhline(0,color='k',linestyle='dashed')

	ax3a.set_xlim(-10,40)
	ax3a.set_ylim(-0.005,0.005)

	# arrows
	xy1 = (ax3a.get_xlim()[0], ax3a.get_ylim()[1])
	xy2 = (-0.2, ax2a.get_ylim()[0])
	col = d_black#(0.1,0.65,0.1)
	alpha = 0.9
	coordsA = 'data'
	coordsB = 'data'
	con = ConnectionPatch(xy1, xy2, coordsA, coordsB,
				arrowstyle="-", shrinkB=0,
				axesA=ax3a, axesB=ax2a, mutation_scale=10,
				ec=col,fc=col,lw=0.75,alpha=alpha)
	ax3a.add_artist(con)

	xy1 = (ax3a.get_xlim()[1], ax3a.get_ylim()[1])
	xy2 = (0.2, ax2a.get_ylim()[0])
	col = d_black#(0.1,0.65,0.1)
	alpha = 1
	coordsA = 'data'
	coordsB = 'data'
	con = ConnectionPatch(xy1, xy2, coordsA, coordsB,
				arrowstyle="-", shrinkB=0,
				axesA=ax3a, axesB=ax2a, mutation_scale=10,
				ec=col,fc=col,lw=0.75,alpha=alpha)
	ax3a.add_artist(con)

	# find zero crossings
	zero_crossing_index = np.where(np.diff(np.signbit(S3_on_zoom.real)))[0]
	zero_crossing_pos_on = d_zoom[zero_crossing_index]
	print('Zero crossing, Thermally populated', zero_crossing_pos_on)
	zero_crossing_index = np.where(np.diff(np.signbit(S3_off_zoom.real)))[0]
	zero_crossing_pos_off = d_zoom[zero_crossing_index]
	print('Zero crossing, Equally populated', zero_crossing_pos_off)

	ax3a.axvline(zero_crossing_pos_on,color=d_grey)
	ax3a.axvline(zero_crossing_pos_off,color=d_grey)

	##
	##
	## argmin (abs (S3)_zoom ) ??



	# fig2.savefig('Boltzmann_onoff.png')
	# fig2.savefig('Boltzmann_onoff.pdf')

	plt.show()

if __name__ == '__main__':
	main()
