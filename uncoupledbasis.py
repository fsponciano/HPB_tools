""" --- Atomic structure calculator in the uncoupled basis ---

	Written by D. J. Whiting (2017)
	** Modified by F. S. Ponciano Ojeda (2021) **
	Quantum Light & Matter Group, Durham University, UK

	Calculates the energy levels of a given state in an alkali-metal
	atom by diagonalising the Hamiltonian in the uncoupled momenta basis,
	suitable for looking at atoms in the Hyperfine Paschen-Back regime.
"""

from numpy import matrix,arange,zeros,kron,eye,pi,dot,linspace,exp,array,sqrt
from scipy.linalg import eigh
from cycler import cycler

import sys,os
# Modify the following to create a directory for saving results of calculations
try:
	os.mkdir('./UncoupledBasis_results')
except OSError:
	pass
try:
	os.mkdir('./UncoupledBasis_figures')
except OSError:
	pass

import matplotlib.pyplot as plt
import fractions
import atomdata
import numpy as np

from scipy.constants import *

h = 2*pi*hbar
muB = physical_constants['Bohr magneton'][0]
kB = physical_constants['Boltzmann constant'][0]
e0 = epsilon_0 #Permittivity of free space
a0 = physical_constants['Bohr radius'][0]
S = 1/2.0 #Electron spin
gS = -physical_constants['electron g factor'][0]

from elecsus.libs.durhamcolours import * # Optional: use colour palette from ElecSus installation

plt.rc('font',**{'family':'Serif','serif':['Times New Roman'],'weight':'bold'})
params={'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14,'legend.fontsize':12,'mathtext.fontset':'cm','mathtext.rm':'serif'}
plt.rcParams.update(params)

def AM_Matrix(l):
	''' Builds a 2l + 1 square matrix for a general angular momentum *l* based
	on construction of the raising and lowering operators. '''

	dim = int(2*l + 1.0)
	m = arange(l,-l - 1,-1)
	LPlus = array(zeros((dim,dim)))
	lm = sqrt((l*(l + 1) - m*(m + 1)))
	LPlus[range(0,dim - 1),range(1,dim)] = lm[1:]
	LMinus = LPlus.transpose()
	Lx = (1/2.0)*(LPlus + LMinus)
	Ly = (-1j/2.0)*(LPlus - LMinus)
	Lz = (1/2.0)*(LPlus @ LMinus - LMinus @ LPlus)
	LL = Lx,Ly,Lz
	return LL

def AM_MatAdd(j1,j2):
	''' Adds two matrices by doing the tensor product of each of the component
	square matrices with their counterpart's identity matrix. '''

	dim1 = j1[0].shape[0]
	dim2 = j2[0].shape[0]
	Jx = kron(j1[0],eye(dim2)) + kron(eye(dim1),j2[0])
	Jy = kron(j1[1],eye(dim2)) + kron(eye(dim1),j2[1])
	Jz = kron(j1[2],eye(dim2)) + kron(eye(dim1),j2[2])
	JJ = Jx,Jy,Jz
	return JJ

def EvalsEvecs(I,L,S,J=0,B=0,IsotopeShift=0,PB=False):
	''' Calculates the eigenvalues (energies) of the atomic system.
	Valid for linear Zeeman and hyperfine Paschen-Back regimes. '''

	SS,LL,II = AM_Matrix(S),AM_Matrix(L),AM_Matrix(I)

	# The order in which angular momenta are combined here defines the order
	# of the decoupled angular momenta in the eigenvectors.
	# Current order gives [mL[mS[mI]]]
	JJ = AM_MatAdd(LL,SS)
	FF = AM_MatAdd(JJ,II)

	S2 = SS[0] @ SS[0] + SS[1] @ SS[1] + SS[2] @ SS[2]
	L2 = LL[0] @ LL[0] + LL[1] @ LL[1] + LL[2] @ LL[2]
	J2 = JJ[0] @ JJ[0] + JJ[1] @ JJ[1] + JJ[2] @ JJ[2]
	I2 = II[0] @ II[0] + II[1] @ II[1] + II[2] @ II[2]
	F2 = FF[0] @ FF[0] + FF[1] @ FF[1] + FF[2] @ FF[2]

	dimS = SS[0].shape[0]
	dimL = LL[0].shape[0]
	dimJ = JJ[0].shape[0]
	dimI = II[0].shape[0]
	dimF = FF[0].shape[0]

	### Set up the hamiltonian in the uncoupled (S,L,I) basis ###
	H = 0
	# Import atomic constants - currently only for Rb
	A_fs, A_hfs, B_hfs, gI, gL = atomdata.atomic_structure_coefficients('Rb',I,L,J)

	# Fine structure
	LdotS = 1/2.0 * kron(J2 - (kron(L2,eye(dimS)) + kron(eye(dimL),S2)),eye(dimI))

	# Spin orbit interaction
	if J == None:
		H += A_fs * LdotS
	else:
		# Recenter the appropriate J state on zero energy
		so_correction_factor = -(J*(J + 1) - L*(L + 1) - S*(S + 1))/2.0
		H += A_fs * (LdotS + so_correction_factor*eye(dimF))

	# Hyperfine interaction
	IdotJ = 1/2.0 * (F2 - (kron(eye(dimJ),I2) + kron(J2,eye(dimI))))

	# Magnetic dipole interaction for hyperfine structure
	H += A_hfs * IdotJ

	# Electric quadrupole interaction for hyperfine structure
	if B_hfs != 0:
		H += B_hfs * (3*(IdotJ @ IdotJ) + 3/2.0*IdotJ - I*(I + 1)*J*(J + 1)*eye(dimF))/(2*I*(2*I - 1)*J*(2*J - 1))

	# Zeeman interaction (assumed along z-axis)
	H -= muB*B/h * (gL*kron(kron(LL[2],eye(dimS)),eye(dimI)) + gS*kron(kron(eye(dimL),SS[2]),eye(dimI)) + gI*kron(kron(eye(dimL),eye(dimS)),II[2]))

	# Isotope shift
	H += IsotopeShift * eye(dimF)

	# Return the eigenvalues and eigenvectors of the hamiltonian (eigensystem)
	return eigh(H)

def BreitRabi(I,L,S,J=None,Bmax=1.5,ylim=None,save=False):
	''' Uses the eigenvalues of the atomic system to generate a Breit-Rabi
	(Zeeman shift) diagram. '''

	B_vals = linspace(0,Bmax,int(3e3))
	degen = (2*I + 1)*(2*L + 1)*(2*S + 1)
	E_vals = zeros((len(B_vals),int(degen)))
	for i in range(len(B_vals)):
		E_vals[i,:] = EvalsEvecs(I,L,S,J,B_vals[i])[0] #Energies, in Hz

	fig = plt.figure(facecolor='None',figsize=(10,8))
	ax = fig.add_subplot(111)

	# Optional use of custom colour-cycle given by m_I value for the ground state
	plot_colours = cycler(color=[d_blue,d_red,d_purple,d_olive,d_midblue,d_pink,d_lightpurple,d_grey,d_black])
	ax.set_prop_cycle(plot_colours[:int(2*I + 1)].concat(plot_colours[:int(2*I + 1)][::-1]))

	ax.plot(B_vals,E_vals/1e9,lw=2)
	ax.set_xlabel('Magnetic field strength (T)',fontweight='bold')
	ax.set_ylabel(r'Energy / $\hbar$ (GHz)',fontweight='bold')
	ax.set_xlim(0,Bmax)
	if ylim != None:
		ax.set_ylim(-ylim - 1,ylim + 1)
	elif ylim == None:
		ax.set_ylim(E_vals[-1].min()/1e9 - 1,E_vals[-1].max()/1e9 + 1)
	fig.subplots_adjust(left=0.09,right=0.95,bottom=0.15,top=0.95)

	# Optional text for adding details of the atom plotted
	fig.text(0.245,0.81,r'$^{87}$Rb',size=16,ha='center',fontweight='bold')
	fig.text(0.245,0.785,r'L = '+str(L)+', J = '+str(fractions.Fraction(J))+', I = '+str(fractions.Fraction(I)),size=16,ha='center',fontweight='bold')

	# Place an indicator at the field at which the HPB regime is entered
	if ylim != None:
		if np.logical_and(L == 0,(2*J + 1) == 2):
			ax.vlines(0.49,-ylim-1,ylim+1,color='r',ls=':')
		if np.logical_and(L == 1,(2*J + 1) == 2):
			ax.vlines(0.06,-ylim-1,ylim+1,color='r',ls=':')
		elif np.logical_and(L == 1,(2*J + 1) == 4):
			ax.vlines(0.04,-ylim-1,ylim+1,color='r',ls=':')
	else:
		if np.logical_and(L == 0,(2*J + 1) == 2):
			ax.vlines(0.49,E_vals[-1].min()/1e9 - 1,E_vals[-1].max()/1e9 + 1,color='r',ls=':')
		if np.logical_and(L == 1,(2*J + 1) == 2):
			ax.vlines(0.06,E_vals[-1].min()/1e9 - 1,E_vals[-1].max()/1e9 + 1,color='r',ls=':')
		elif np.logical_and(L == 1,(2*J + 1) == 4):
			ax.vlines(0.04,E_vals[-1].min()/1e9 - 1,E_vals[-1].max()/1e9 + 1,color='r',ls=':')

	fig.text(0.135,0.675,r'$(F,m_{F})$',size=16,ha='center',fontweight='bold')
	fig.text(0.865,0.705,r'$(J,m_{J};I,m_{I})$',size=16,ha='center',fontweight='bold')

	plt.show()

	if save == True:
		plt.savefig('./UncoupledBasis_figures/Bmax_'+str(B)+'_Rb_I-'+str(I)+'_L-'+str(L)+'_J-'+str(J)+'.pdf',dpi=300)
		plt.savefig('./UncoupledBasis_figures/Bmax_'+str(B)+'_Rb_I-'+str(I)+'_L-'+str(L)+'_J-'+str(J)+'.png',dpi=300)

def AM_StateDecomp(I,L,S,J=None,B=0,ylim=None,cutoff=0.001):
	''' Calculates the admixture/decomposition of the atomic energy levels
	in terms of the uncoupled basis (m_L,m_S,m_I) states. '''

	output_states = []
	if ylim == None:
		ylim = np.inf
	Atom_eigsys = EvalsEvecs(I,L,S,J,B)
	state_eigvecs = abs(Atom_eigsys[1])
	state_energies = Atom_eigsys[0]/1e9
	state_eigvecs,state_energies = state_eigvecs[:,::-1],state_energies[::-1] # Reverse the array so highest energy comes first
	print('Decoupled quantum numbers with fractions >', cutoff)
	for i in range(state_eigvecs.shape[1]):
		if abs(state_energies[i]) < ylim + 1:
			print('State ',(i+1))
			evec = state_eigvecs[:,i]
			print('----------------------')
			print('Energy / hbar (GHz):', state_energies[i])
			print('Fractional composition of |mI,  mL,  mS > states:')
			for iI,mI in enumerate(arange(-I,I + 1)):
				for iL,mL in enumerate(arange(-L,L + 1)):
					for iS,mS in enumerate(arange(-S,S + 1)):
						frac = evec[int(iL*(2*I + 1)*(2*S + 1) + iS*(2*I + 1) + iI)]
						if frac > cutoff:
							state_decomp = '%.4f' % state_energies[i], '%.4f' % frac, r'$| %2s, %2s, %4s >$' % (str(fractions.Fraction(mI)),str(fractions.Fraction(mL)),str(fractions.Fraction(mS)))
							print(state_decomp)
							output_states.append(state_decomp)
							#print mL+mI+mS
			print('----------------------\n')
			# np.savetxt('./Uncoupled basis results/DecoupledStates_'+str(I)+'_'+str(L)+'_'+str(J)+'_'+str(B)+'T.csv',output_states,delimiter=',')
	return output_states

if __name__=="__main__":
	I = 3/2
	L = 0
	S = 1/2
	J = 1/2
	Bvalues = [0.2,0.6,1.5] #in Tesla

	print('Calculating energies for B =',Bvalues)
	for B in Bvalues:
		print('B = '+str(B))
		print('Energies: \n')
		print(EvalsEvecs(I,L,S,J,B)[0])
		print('-----\t-----\t-----')
		## Uncomment line below to save raw eigenvectors/values for the system
		# np.savetxt('./UncoupledBasis_results/EigenVals_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'_'+str(B)+'T.csv',EvalsEvecs(I,L,S,J,B),delimiter=',')

		g = 3/2.0 + (S*(S + 1) - L*(L + 1))/(2*J*(J + 1))
		ylim = 1.5*J*((g*muB*B/h)/1e9)
		BreitRabi(I,L,S,J,B,ylim=ylim,save=False)
		decomp_states = AM_StateDecomp(I,L,S,J,B,ylim=ylim)
		## Uncomment lines below to save state decompositions
		# results_decomp = './UncoupledBasis_results/StateDecomp_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'_'+str(B)+'T.txt'
		# with open(results_decomp,'w') as file_output:
		# 	print >> file_output, AM_StateDecomp(I,L,S,J,B,ylim=ylim)

	plt.show()
