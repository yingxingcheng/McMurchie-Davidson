import time
import numpy as np

np.set_printoptions(precision=6, linewidth=np.inf, suppress=True)
import psi4

# Memory for Psi4 in GB
psi4.set_memory("500 MB")
psi4.core.set_output_file("output.dat", False)

# Memory for numpy in GB
numpy_memory = 2

mol = psi4.geometry(
    """
O 0 0 0.0
O 0 0 0.74
symmetry c1
units angstrom
no_com
no_reorient
"""
)

# bs_fname = "car_sto-3g"
# bs_fname = "my_sto-3g"
# bs_fname = "psi4-3-21g"
# bs_fname = "psi4-sto-3g"
bs_fname = "psi4-cc-pvdz"
# bs_fname = "sto-3g"

# psi4.set_options({"basis": "cc-pvdz", "scf_type": "pk", "e_convergence": 1e-8})
# psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "e_convergence": 1e-8})
# psi4.set_options({"basis": "my_sto-3g", "scf_type": "pk", "e_convergence": 1e-8})
psi4.set_options({"basis": bs_fname, "scf_type": "pk", "e_convergence": 1e-8})

# Integral generation from Psi4's MintsHelper
wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
# wfn = psi4.core.Wavefunction.build(mol, bs_fname)
t = time.time()
mints = psi4.core.MintsHelper(wfn.basisset())
S = np.asarray(mints.ao_overlap())

# Get nbf and ndocc for closed shell molecules
nbf = S.shape[0]
# ndocc = wfn.nalpha()

# print("\nNumber of occupied orbitals: %d" % ndocc)
# print("Number of basis functions: %d" % nbf)

# # Compute required quantities for SCF
# V = np.asarray(mints.ao_potential())
# T = np.asarray(mints.ao_kinetic())
# # I = np.asarray(mints.ao_eri())

# S_horton = np.array(
#     [
#         [1.0, 0.236704, 0.0, 0.0, 0.0, 0.000381, 0.094791, 0.0, 0.0, 0.160385],
#         [0.236704, 1.0, 0.0, 0.0, 0.0, 0.094791, 0.609404, 0.0, 0.0, 0.505654],
#         [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.437992, 0.0, 0.0],
#         [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.437992, 0.0],
#         [0.0, 0.0, 0.0, 0.0, 1.0, -0.160385, -0.505654, 0.0, 0.0, -0.196673],
#         [0.000381, 0.094791, 0.0, 0.0, -0.160385, 1.0, 0.236704, 0.0, 0.0, 0.0],
#         [0.094791, 0.609404, 0.0, 0.0, -0.505654, 0.236704, 1.0, 0.0, 0.0, 0.0],
#         [0.0, 0.0, 0.437992, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
#         [0.0, 0.0, 0.0, 0.437992, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
#         [0.160385, 0.505654, 0.0, 0.0, -0.196673, 0.0, 0.0, 0.0, 0.0, 1.0],
#     ]
# )

print(S)
print(S.shape)
# print(S / S_horton)

np.savez(
    "psi4_oxygen_pairs_overlap_2d_cc_pvdz.npz",
    xyz=np.asarray(mol.geometry()),
    overlap=S,
)
