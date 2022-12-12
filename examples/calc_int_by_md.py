from mmd.molecule import Molecule
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

np.set_printoptions(suppress=True, precision=6, linewidth=np.inf)

# Simple example of Born-Oppenheimer Molecular Dynamics
# using minimal basis H2. Shows how you can chain together
# routines to create more complex programs. This runs about
# a a half of a femtosecond of molecular dynamics and plots
# energy as a function of time

# Molecular geometry input
h2 = """
0 1
O 0.0 0.0 0.00
O 0.0 0.0 0.74
"""

# h2 = """
# 0 1
# Li 0.0 0.0 0.00
# Li 0.0 0.0 0.7
# """

# init molecule and build integrals
# mol = Molecule(geometry=h2, basis="sto-3g")
# mol = Molecule(geometry=h2, basis="gauss-sto-3g")
mol = Molecule(geometry=h2, basis="gauss-cc-pvdz")
mol.build()
S = mol.S
print(S)
print(S.shape)

# trans = np.array(
#     [
#         [-1 / 2, 0, 0, -1 / 2, 0, 1],
#         [0, 0, 1, 0, 0, 0],
#         [0, 0, 0, 0, 1, 0],
#         [1 / 2 * np.sqrt(3), 0, 0, -1 / 2 * np.sqrt(3), 0, 0],
#         [0, 1, 0, 0, 0, 0],
#     ]
# )
# print(trans)
#
# tf2 = np.array(
#     [
#         [-0.5, 0, 0, -0.5, 0, 1.0],
#         [0, 0, 1.7320508075688773, 0, 0, 0],
#         [0, 0, 0, 0, 1.7320508075688773, 0],
#         [0.86602540378443865, 0, 0, -0.86602540378443865, 0, 0],
#         [0, 1.7320508075688773, 0, 0, 0, 0],
#     ]
# )
# print(tf2)
#
# print("overlap with pure Gaussian")
# print(tf2 @ S @ tf2.T)
# print(trans @ S @ trans.T)
