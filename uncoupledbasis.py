""" --- Atomic structure calculator in the uncoupled basis --- """
""" Written by D. J. Whiting (2017);
	** Modified by F. S. Ponciano Ojeda (2021) **
	@ QLM, Durham University """
""" Calculates the energy levels of an atomic transition in an alkali-metal
atom by diagonalising the Hamiltonian in the uncoupled momenta basis, suitable
for looking at the Hyperfine Paschen-Back regime. """

# from __future__ import division,print_function # Uncomment for Python 2 compatibility
from numpy import matrix,arange,zeros,kron,eye,pi,dot,linspace,exp,array
from scipy.linalg import eigh

import sys
sys.path.append("C:\'Users\'Francisco\'Documents\'Voigt_Boltzmann\'kB_Big_Diagram") #Change for general path

import matplotlib.pyplot as plt
import matplotlib
import fractions
import atomdata
import numpy as np

from scipy.constants import *
h = 2*pi*hbar
muB = physical_constants['Bohr magneton'][0]
kB = physical_constants['Boltzmann constant'][0]
e0 = epsilon_0 #Permittivity of free space
a0 = physical_constants['Bohr radius'][0]
S = 0.5 #Electron spin
gS = -physical_constants['electron g factor'][0]

from durhamcolours import *

plt.rc('font',**{'family':'Serif','serif':['Times New Roman'],'weight':'bold'})
params={'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14,'legend.fontsize':12,'mathtext.fontset':'cm','mathtext.rm':'serif'}
plt.rcParams.update(params)

def jmat(j):
	''' Builds a 2j+1 square matrix for the total orbital angular momentum,
	defined as J = L + S. '''

	dim = int(2*j + 1.0)
	m = arange(j,-j - 1,-1)
	JPlus = matrix(zeros((dim,dim)))
	jm = (j*(j + 1) - m*(m + 1))**(1/2.0)
	JPlus[range(0,dim - 1),range(1,dim)] = jm[1:]
	JMinus = JPlus.transpose()
	Jx = (1/2.0)*(JPlus + JMinus)
	Jy = (-1j/2.0)*(JPlus - JMinus)
	Jz = (1/2.0)*(JPlus*JMinus - JMinus*JPlus)
	JJ = Jx,Jy,Jz
	return JJ

def jsum(j1,j2):
	''' Adds matrices by doing the tensor product of each of the component
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

	SS,LL,II = jmat(S),jmat(L),jmat(I)

	### The order in which angular momenta are combined here defines the order
	### of the decoupled angular momenta in the eigenvectors.
	### Current order gives [mL[mS[mI]]]
	JJ = jsum(LL,SS)
	FF = jsum(JJ,II)

	S2 = SS[0]**2 + SS[1]**2 + SS[2]**2
	L2 = LL[0]**2 + LL[1]**2 + LL[2]**2
	J2 = JJ[0]**2 + JJ[1]**2 + JJ[2]**2
	I2 = II[0]**2 + II[1]**2 + II[2]**2
	F2 = FF[0]**2 + FF[1]**2 + FF[2]**2

	dimS = SS[0].shape[0]
	dimL = LL[0].shape[0]
	dimJ = JJ[0].shape[0]
	dimI = II[0].shape[0]
	dimF = FF[0].shape[0]

	#### Set up the hamiltonian in the uncoupled (S,L,I) basis ####
	H = 0
	# Import atomic constants - currently only for Rb
	A_fs, A_hfs, B_hfs, gI, gL = atomdata.atomic_structure_coefficients('Rb',I,L,J)

	# Fine structure
	LdotS = 1/2 * kron(J2 - (kron(L2,eye(dimS)) + kron(eye(dimL),S2)),eye(dimI))

	# Spin orbit interaction
	if J == None:
		H += A_fs * LdotS
	else:
		# Recenter the appropriate J state on zero energy
		so_correction_factor = -(J*(J + 1) - L*(L + 1) - S*(S + 1))/2
		H += A_fs * (LdotS + so_correction_factor*eye(dimF))

	# Hyperfine structure
	IdotJ = 1/2 * (F2 - (kron(eye(dimJ),I2) + kron(J2,eye(dimI))))

	# Magnetic dipole interaction
	H += A_hfs * IdotJ

	# Electric quadrupole interaction
	if B_hfs != 0:
		H += B_hfs * (3*IdotJ**2 + 3/2*IdotJ - I*(I + 1)*J*(J + 1)*eye(dimF))/(2*I*(2*I - 1)*J*(2*J - 1))

	# Zeeman interaction
	H -= muB*B/h * (gL*kron(kron(LL[2],eye(dimS)),eye(dimI)) + gS*kron(kron(eye(dimL),SS[2]),eye(dimI)) + gI*kron(kron(eye(dimL),eye(dimS)),II[2]))

	# Isotope shift
	H += IsotopeShift * eye(dimF)

	return eigh(H) # Return the eigenvalues and eigenvectors of the hamiltonian

def BreitRabi(I,L,S,J=None,Bmax=1.5,ylim=None,save=False):
	''' Uses the eigenvalues of the atomic system to generate a Breit-Rabi
	(Zeeman shift) diagram. '''

	Bs = linspace(0,Bmax,int(3e3))
	d = (2*I + 1)*(2*L + 1)*(2*S + 1)
	Es = zeros((len(Bs),int(d)))
	for i in range(len(Bs)):
		Es[i,:] = EvalsEvecs(I,L,S,J,Bs[i])[0]

	fig = plt.figure(facecolor='white', figsize=(10,8))
	ax = fig.add_subplot(111)
	ax.plot(Bs,Es/1e9,'k')
	ax.set_xlabel('Magnetic field strength (T)',fontweight='bold')
	ax.set_ylabel('Energy (GHz)',fontweight='bold')
	ax.set_xlim(0,Bmax)
	if ylim != None:
		ax.set_ylim(-ylim-1,ylim+1)
	# fig.text(0.25,0.775,r'B = '+str(Bmax)+' T',size=16,ha='center')
	fig.subplots_adjust(left=0.09,right=0.95,bottom=0.15,top=0.95)

	fig.text(0.25,0.8,r'$^{87}$Rb',size=16,ha='center',fontweight='bold')
	fig.text(0.25,0.775,r'L = '+str(L)+', J = '+str(fractions.Fraction(J))+', I = '+str(fractions.Fraction(I)),size=16,ha='center',fontweight='bold')

	if np.logical_and(L == 0,(2*J+1) == 2):
		ax.vlines(0.49,-ylim-1,ylim+1,color='r',ls=':')
	if np.logical_and(L == 1,(2*J+1) == 2):
		ax.vlines(0.06,-ylim-1,ylim+1,color='r',ls=':')
	elif np.logical_and(L == 1,(2*J+1) == 4):
		ax.vlines(0.04,-ylim-1,ylim+1,color='r',ls=':')

	fig.text(0.135,0.675,r'$(F,m_{F})$',size=16,ha='center',fontweight='bold')
	fig.text(0.865,0.895,r'$(J,m_{J};I,m_{I})$',size=16,ha='center',fontweight='bold')

	if save == True:
		plt.savefig('./Plots/B_'+str(B)+'_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'.pdf',dpi=300)
		plt.savefig('./Plots/B_'+str(B)+'_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'.png',dpi=300)
		# np.savetxt('./BR_Eigen_'+str(B)+'_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'.csv',Es[0,:],delimiter=',')

def MomentumDecomp(I,L,S,J=None,B=0,ylim=None,cutoff=0.001):
	''' Calculates the admixture/decomposition of the atomic energy levels
	in terms of the uncoupled basis (m_L,m_S,m_I) states. '''
	output_states = []
	if ylim == None:
		ylim = np.inf
	evalsevecs = EvalsEvecs(I,L,S,J,B)
	evecs = abs(evalsevecs[1])
	evals = evalsevecs[0]/1e9
	evecs,evals=evecs[:,::-1],evals[::-1] # Reverse the array so highest energy comes first
	print('Decoupled quantum numbers with fractions >', cutoff)
	for i in range(evecs.shape[1]):
		if abs(evals[i])<ylim+1:
			evec = evecs[:,i]
			print('----------------------')
			print('Energy:', evals[i])
			print('Fraction		|mI,  mL,  mS >')
			for iI,mI in enumerate(arange(-I,I+1)):
				for iL,mL in enumerate(arange(-L,L+1)):
					for iS,mS in enumerate(arange(-S,S+1)):
						frac = evec[int(iL*(2*I+1)*(2*S+1)+iS*(2*I+1)+iI)]
						if frac > cutoff:
							state_decomp = '%.4f' % evals[i], '%.4f' % frac, r'$| %2s, %2s, %4s >$' % (str(fractions.Fraction(mI)),str(fractions.Fraction(mL)),str(fractions.Fraction(mS)))
							print(state_decomp)

							output_states.append(state_decomp)
							#print mL+mI+mS
			print('----------------------')
			# np.savetxt('./Uncoupled basis results/DecoupledStates_'+str(I)+'_'+str(L)+'_'+str(J)+'_'+str(B)+'T.csv',output_states,delimiter=',')
	return output_states

if __name__=="__main__":
	I = 3/2
	L = 1
	S = 1/2
	J = 3/2
	# Bvalues = [0.4,0.6,1.5,8,15,35] # Tesla
	Bvalues = [0.2,0.6,1.5]

	for B in Bvalues:
		results_decomp = './Uncoupled basis results/StateDecomp_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'_'+str(B)+'T.txt'
		print('B = '+str(B))
		print('Energies: \n')
		print(EvalsEvecs(I,L,S,J,B)[0])
		print('-----\t-----\t-----')
		# np.savetxt('./Uncoupled basis results/EigenVals_Rb_'+str(I)+'_'+str(L)+'_'+str(J)+'_'+str(B)+'T.csv',EvalsEvecs(I,L,S,J,B)[0],delimiter=',')

		g = 3/2 + (S*(S+1)-L*(L+1))/(2*J*(J+1))
		ylim = 1.5*J*((g*muB*B/h)/1e9)
		BreitRabi(I,L,S,J,B,ylim=ylim,save=False)
		MomentumDecomp(I,L,S,J,B,ylim=ylim)
		# with open(results_decomp,'w') as file_output:
		# 	print >> file_output, MomentumDecomp(I,L,S,J,B,ylim=ylim)

	plt.show()
