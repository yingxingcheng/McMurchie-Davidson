import numpy as np
from horton import *

log.set_level(0)

np.set_printoptions(suppress=True, precision=6, linewidth=np.inf)

mol = IOData.from_file("mol.xyz")

# obasis = get_gobasis(mol.coordinates, mol.numbers, "sto-3g", element_map={1: "3-21g", 6: "6-31g"})
# obasis = get_gobasis(mol.coordinates, mol.numbers, "sto-3g")
# obasis = get_gobasis(mol.coordinates, mol.numbers, "6-31gss")
# obasis = get_gobasis(mol.coordinates, mol.numbers, "6-31g")

mybasis = GOBasisFamily("my-cc-pvdz", filename="horton-cc-pvdz.nwchem")
is_pure = True
obasis = get_gobasis(mol.coordinates, mol.numbers, mybasis, pure=is_pure)

# Create a linalg factory
lf = DenseLinalgFactory(obasis.nbasis)

# Compute Gaussian integrals
olp = obasis.compute_overlap(lf)
# kin = obasis.compute_kinetic(lf)
# na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
# er = obasis.compute_electron_repulsion(lf)


print(obasis.nbasis)
print(olp._array)
# print(kin._array)
print(olp._array.shape)

if not is_pure:
    np.savez(
        "horton_oxygen_pairs_overlap_2d_cc_pvdz_cart.npz",
        xyz=mol.coordinates,
        numbers=mol.numbers,
        overlap=olp._array,
    )
else:
    np.savez(
        "horton_oxygen_pairs_overlap_2d_cc_pvdz.npz",
        xyz=mol.coordinates,
        numbers=mol.numbers,
        overlap=olp._array,
    )
