""" --- Atom data ---

	Written by D. J. Whiting (2017)
	** Modified by F. S. Ponciano Ojeda (2021) **
	Quantum Light & Matter Group, Durham University, UK

    Simple database for housing atomic species data for use in calculations of
    energy levels.
"""
def atomic_structure_coefficients(atom,I,L,J):
    '''
        Currently only contains information for Rb states with n=5, l=0,1,2
        and Cs states with n=6, l=0,1,2. Further data to be added...
    '''
    # global  A_fs, A_hfs, B_hfs, gI, gL
    if atom == 'Rb':
        if int(2*I+1e-10) == 3: # Rb 87
            gI = -0.0009951414
            gL = 0.99999369
            if L == 0:
                A_fs = 0
                A_hfs = 3417.34130545215e6
                B_hfs = 0
            elif L == 1:
                A_fs = 7.123e12*2/3 # (2*S+1)/(2*L+1) or 1/(S+L)
                if int(2*J+1e-10) == 1:
                    A_hfs = 406.147e6
                    B_hfs = 0
                elif int(2*J+1e-10) == 3:
                    A_hfs = 84.7185e6
                    B_hfs = 12.4965e6
            elif L == 2:
                A_fs = 88.905185e9*2/5 # (2*S+1)/(2*L+1) or 1/(S+L)
                if int(2*J+1e-10) == 3:
                    A_hfs = 14.4303e6
                    B_hfs = 0.9320e6
                elif int(2*J+1e-10) == 5:
                    A_hfs = -7.4605e6
                    B_hfs = 1.2713e6
        elif int(2*I+1e-10) == 5: # Rb 85
            gI = -0.0002936400
            gL = 0.99999354
            if L == 0:
                A_fs = 0
                A_hfs = 1011.9108130e6
                B_hfs = 0
            elif L == 1:
                A_fs = 7.123e12*2/3
                if int(2*J+1e-10) == 1:
                    A_hfs = 120.527e6
                    B_hfs = 0
                elif int(2*J+1e-10) == 3:
                    A_hfs = 25.0020e6
                    B_hfs = 25.790e6
            elif L == 2:
                A_fs = 88.905185e9*2/5
                if int(2*J+1e-10) == 3:
                    A_hfs = 4.2699e6
                    B_hfs = 1.9106e6
                elif int(2*J+1e-10) == 5:
                    A_hfs = -2.2112e6
                    B_hfs = 2.6804e6

    if atom == 'Cs':
        gI = -0.00039885395
        gL = 0.99999587
        if L == 0:
            A_fs = 0
            A_hfs = 2298.1579425e6
            B_hfs = 0
        elif L == 1:
            A_fs = 16.60966969094428436e12*2/3 # (2*S+1)/(2*L+1) or 1/(S+L)
            if int(2*J+1e-10) == 1:
                A_hfs = 291.9201e6
                B_hfs = 0
            elif int(2*J+1e-10) == 3:
                A_hfs = 50.28827e6
                B_hfs = −0.4934e6
                C_hfs = 0.560e3
        elif L == 2:
            A_fs = 1.285060371217e12*2/5 # (2*S+1)/(2*L+1) or 1/(S+L)
            if int(2*J+1e-10) == 3:
                A_hfs = 16.34e6
                B_hfs = -0.1e6
            elif int(2*J+1e-10) == 5:
                A_hfs = -4.66e6
                B_hfs = 0.9e6

    return A_fs, A_hfs, B_hfs, gI, gL
