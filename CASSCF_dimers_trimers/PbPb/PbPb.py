from pyscf import gto, scf, ao2mo, fci, mcscf, tools, cc
from pyscf.lib import chkfile
import numpy as np

bond_lengths = np.linspace(2.2, 14.0,119)
a = list(bond_lengths)
#a.reverse()
print(bond_lengths)
k = 0
#energy = []
f = open('track_CAS_PbCs_1010.dat', 'w')
for bond_length in a: 
    print("Bond-length: "+str(bond_length))
    mol = gto.Mole()
    mol.atom = [['Pb',(0,0,0)],['Pb',(0,0,bond_length)]]
    mol.basis = 'def2-tzvpp'
    mol.ecp = 'def2-tzvpp'
    mol.symmetry = 'D2h' 
    mol.build()
    mf = scf.RHF(mol)
    mf.chkfile = 'Pbchk'
    
    mf.verbose = 3
    mf.max_cycle = 10000
    #mf.conv_tol = 1e-6
    if k == 0:
        mf.kernel()
    #    k = 1
    else:
        mf.init_guess = 'chkfile'
        mf.kernel()
    mc = mcscf.CASSCF(mf, 10, 10)
    if k == 0:
        cas_e = mc.kernel()[0]
        cas_mo = mc.mo_coeff.copy()
        cas_ci = mc.ci.copy()
        k = 1
    else:
        mo_init = mcscf.project_init_guess(mc, cas_mo)
        cas_e = mc.kernel(mo_init, cas_ci)[0]
        cas_ci = mc.ci.copy()
        cas_mo = mc.mo_coeff.copy()
#    mcc = cc.CCSD(mf).run()
#    et = mcc.ccsd_t()
    f.write(str(bond_length) + '    ' + str(mf.e_tot) + '    ' + str(cas_e) + '\n')
f.close()
    
#    print(mf.mo_coeff)
#    tools.fcidump.from_mo(mol, 'FCIDUMP.PbI2.def2-tzvpp-'+str(bond_length), mf.mo_coeff)
    #tools.fcidump.from_mo(mol, 'CASSCF-dump', mc.mo_coeff)
    
    #h1 = mf.mo_coeff.T.dot(mf.get_hcore()).dot(mf.mo_coeff)
    #eri = ao2mo.kernel(mol, mf.mo_coeff)
    #cisolver = fci.direct_spin1.FCI(mol)
    #e, ci = cisolver.kernel(h1,eri,h1.shape[1], mol.nelec, ecore= mol.energy_nuc())
    #print(e)
