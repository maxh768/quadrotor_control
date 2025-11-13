import numpy as np
from dynamics import *
from scipy.linalg import solve_continuous_are

def lqr(A, B, Q, R):
    # Solve the continuous time algebraic Riccati equation
    P = solve_continuous_are(A, B, Q, R)

    # Compute the LQR gain
    K = np.linalg.inv(R) @ B.T @ P

    return K

def form_closed_loop_system(params, x_eq, u_eq, Q, R):
    """
    Compute the closed-loop system matrices A_cl, B_cl
    for the linearized system with LQR controller
    """
    A, B = linearize_system(0, x_eq, u_eq, params)
    K = lqr(A, B, Q, R)

    A_cl = A - B @ K

    return A_cl, K



if __name__ == "__main__":
    m = 1.4
    I = 0.0211
    l = 0.159
    beta = 0.1365
    g = 9.81
    params = (m, I, l, beta)

    x_eq = np.array([0., 0., 0., 0., 0., 0.])
    u_eq = np.array([m * g / 2, m * g / 2])

    Q = np.eye(6)
    R = np.eye(2)

    A_cl = form_closed_loop_system(params, x_eq, u_eq, Q, R)

    print(A_cl)

    A, B = linearize_system(0, x_eq, u_eq, params)
    eig_OL , _ = np.linalg.eig(A)
    print("Open-loop eigenvalues:", eig_OL)

    eig, _ = np.linalg.eig(A_cl)
    print("Closed-loop eigenvalues:", eig)