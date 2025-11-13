import numpy as np

def nonlinear_dynamics(t, state, u, params):
    """
    Full nonlinear dynamics of 2D quadrotor
    """
    # state
    x, y, theta, v_x, v_y, omega = state
    
    # params
    m, I, l, beta = params
    g = 9.81

    # inputs
    T1, T2 = u

    dx = v_x
    dy = v_y
    dtheta = omega
    dv_x = ( -2 * (T1 + T2) * np.sin(theta) \
            - beta * v_x * np.sqrt(v_x**2 + v_y**2) ) / m
    
    dv_y = ( ( 2 * (T1 + T2) * np.cos(theta) \
            - beta * v_y * np.sqrt(v_x**2 + v_y**2) ) / m ) - g
        
    omega_dot = 2 * (T1 - T2) * l / I

    dx = np.array([dx, dy, dtheta, dv_x, dv_y, omega_dot])

    return dx

def linearize_system(t, x, u, params):
    """
    Linearize the nonlinear dynamics around given x, u
    """
    n = x.shape[0]
    m = u.shape[0]

    A = np.zeros((n, n))
    B = np.zeros((n, m))
    eps = 1e-6

    dx = nonlinear_dynamics(0, x, u, params)

    for i in range(n):
        x_p = x.copy()
        x_p[i] += eps

        dx_p = nonlinear_dynamics(0, x_p, u, params)

        der = (dx_p - dx) / eps
        A[:, i] = der

    for i in range(m):
        u_p = u.copy()
        u_p[i] += eps
        dx_p = nonlinear_dynamics(0, x, u_p, params)

        der = (dx_p - dx) / eps
        B[:, i] = der

    return A, B
    
def linear_dynamics(t, x, u, A, B, xbar, ubar):
    """
    Linear dynamics of 2D quadrotor
    """
    dx = A @ (x - xbar) + B @ (u - ubar)
    return dx

def LQR_nonlinear_dynamics(t, x, k, params, xbar, ubar):
    u = - k @ (x - xbar) + ubar
    dx = nonlinear_dynamics(t, x, u, params)
    return dx

if __name__ == "__main__":
    m = 1.4
    I = 0.0211
    l = 0.159
    beta = 0.1365
    g = 9.81
    params = (m, I, l, beta)

    x = np.array([0., 0., 0., 0., 0., 0.])
    u = np.array([m * g / 2, m * g / 2])

    A, B = linearize_system(0, x, u, params)
    print("A:", A)
    print("B:", B)




