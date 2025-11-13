import numpy as np
from dynamics import *
from controller import get_gains


# parameters    
m = 1.4
I = 0.0211
l = 0.159
beta = 0.1365
g = 9.81
params = (m, I, l, beta)

# operating point
x_eq = np.array([0., 0., 0., 0., 0., 0.])
u_eq = np.array([m * g / 2, m * g / 2])

C = np.array([[1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0]])
Q = np.eye(6)
R = np.eye(2) * 0.1

A, B = linearize_system(0, x_eq, u_eq, params)
A_cl, B_cl, K, G = get_gains(params, x_eq, u_eq, Q, R, C)

# stability
eig, _ = np.linalg.eig(A)
max_real_eig = np.max(np.real(eig))
eig_cl, _ = np.linalg.eig(A_cl)
max_real_eig_cl = np.max(np.real(eig_cl))

print(f'Max Real Eigenvalue of A: {max_real_eig}, Max Real Eigenvalue of A_cl: {max_real_eig_cl}')

# controllablity
P_blocks = []
P_blocks.append(B)
for i in range(5):
    A_power = np.linalg.matrix_power(A, i)
    P_block = A_power @ B
    P_blocks.append(P_block)

P = np.hstack(P_blocks)

rank_P = np.linalg.matrix_rank(P)
min_dim = min(P.shape[0], P.shape[1])

if rank_P == min_dim:
    print(f'System is controllable with rank {rank_P}')
else:
    print(f'System is NOT controllable with rank {rank_P}')

# observability
O_blocks = []
O_blocks.append(C)
for i in range(5):
    A_power = np.linalg.matrix_power(A, i)
    O_block = C @ A_power
    O_blocks.append(O_block)

O = np.vstack(O_blocks)

rank_O = np.linalg.matrix_rank(O)
min_dim = min(O.shape[0], O.shape[1])

if rank_O == min_dim:
    print(f'System is observable with rank {rank_O}')
else:
    print(f'System is NOT observable with rank {rank_O}')




