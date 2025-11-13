import numpy as np
from dynamics import *
from scipy.linalg import solve_continuous_are
from scipy.signal import place_poles

def lqr(A, B, Q, R):
    # Solve the continuous time algebraic Riccati equation
    P = solve_continuous_are(A, B, Q, R)

    # Compute the LQR gain
    K = np.linalg.inv(R) @ B.T @ P

    return K

def tracking_gain(A, B, K, C):
    den = C @ np.linalg.inv(B @ K - A) @ B
    G = np.linalg.inv(den)
    return G

def get_gains(params, x_eq, u_eq, Q, R, C):
    """
    Compute the closed-loop system matrices A_cl, B_cl
    for the linearized system with LQR controller
    """
    A, B = linearize_system(0, x_eq, u_eq, params)
    K = lqr(A, B, Q, R)
    G = tracking_gain(A, B, K, C)

    # observer gains
    # need six poles, for now just use purely real poles
    poles_real_part = np.linspace(-10, -15, 6)
    # poles_im_part = 
    L_transpose = place_poles(A.T, C.T, poles_real_part).gain_matrix
    L = L_transpose.T




    A_cl = A - B @ K
    B_cl = B @ G

    return A_cl, B_cl, K, G, L


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

    C = np.array([[1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0]])

    A_cl, B_cl, K, G, L = get_gains(params, x_eq, u_eq, Q, R, C)

    A, B = linearize_system(0, x_eq, u_eq, params)
    eig_OL , _ = np.linalg.eig(A)
    print("Open-loop eigenvalues:", eig_OL)

    eig, _ = np.linalg.eig(A_cl)
    print("Closed-loop eigenvalues:", eig)