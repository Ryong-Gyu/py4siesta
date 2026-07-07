import os, glob
import numpy as np
import matplotlib.pyplot as plt
import math


def delta(x):

    if abs(x) > 8*smearing:
        return 0
    else:
        result = np.exp(-(x/smearing)**2)/(smearing*np.sqrt(pi))
        return result

def readEIG(system_label):


    f = open('%s'%(system_label+'.EIG'))
    tlines = f.readlines()

    for i in range(len(tlines)):
        line = tlines[i]
        words = line.split()
        if i == 0: E_f = float(words[0])     # fermi level
        elif i == 1:
            nbands = int(words[0])           # number of bands
            nspin = int(words[1])            # number of spins
            nk = int(words[2])               # number of k

            k = []
            E = np.zeros((nbands,nk))
        else:
            klabel = (i-2)//int(math.ceil(nbands/10))

            if (i-2)%int(math.ceil(nbands/10))==0:
                k.append(float(words[0]))
                for j in range(1,len(words)):
                    E[j-1][klabel]=float(words[j])

            else:
                blabel = ((i-2)%int(math.ceil(nbands/10)))*10
                for j in range(len(words)):
                    E[blabel+j][klabel]=float(words[j])

    

    VBMs = []
    CBMs = []

    for j in range(nk):
        for i in range(nbands):
            if E[i][j] > E_f:
                VBMs.append(E[i-1][j])
                CBMs.append(E[i][j])
                break

    vbm = max(VBMs)
    cbm = min(CBMs)

    return E, vbm, cbm

def readKP(system_label):

    f = open(system_label + '.KP', 'r')
    lines = f.readlines()
    nk    = int(lines[0].split()[0])

    k     = np.zeros((nk,3), dtype = float)
    wk    = np.zeros(nk, dtype = float)

    for ik in range(nk):
        index = 1+ik
        line = lines[index].split()
        k[ik,0] = float(line[0])
        k[ik,1] = float(line[1])
        k[ik,2] = float(line[2])
        wk[ik] =  float(line[4])

    return k, wk



def dos(nenergy, emin, emax, system_label):

    energy = np.linspace(emin,emax,nenergy)
    eigval, vbm, cbm = readEIG(system_label)
    k, wk = readKP(system_label)
    nbands, nk = eigval.shape
    dos = np.zeros(nenergy, dtype = float)

    for ik in range(nk):
        for ieig in range(nbands):
            eigenvalue = eigval[ieig,ik]
            kweight    = wk[ik]
            for ie in range(nenergy):
                dos[ie] = dos[ie] + kweight * delta(energy[ie]-eigenvalue)

    sumwk = 0
    for ik in range(nk):
        sumwk += wk[ik]

    dos = dos/sumwk

    return energy, dos, vbm, cbm

def cut_range(emin, emax, ef, e, dos):

    ne = len(e)

    energy = e - ef

    DOS = []
    E   = []
    for ie in range(ne):
        if (energy[ie] < emax and emin < energy[ie]):
            DOS.append(dos[ie])
            E.append(energy[ie])

    DOS = np.array(DOS, dtype = float)
    E   = np.array(E, dtype = float)

    return E, DOS

system  = 'Molecular'
emin = -5
emax = 10
smearing = 0.05
pi = np.pi

energy, dos, Ev, Ec  = dos(2000,-15,5,system)
plt.plot(energy, dos, c = 'k')
plt.axvline(x = Ev, c = 'r', alpha = 0.5)
plt.axvline(x = Ec, c = 'b', alpha = 0.5)
plt.xlim([-10,0])

#plt.show()
#E, DOS = cut_range(emin, emax, Ev, energy, dos)
#plt.plot(E,DOS)
plt.savefig('dos.png')
