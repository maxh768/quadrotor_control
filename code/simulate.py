import numpy as np
from dynamics import *
from controller import get_gains
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import niceplots

plt.style.use(niceplots.get_style('james-light'))

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

A, B = linearize_system(0, x_eq, u_eq, params)

# optimal control weighting
Q = np.eye(6)
R = np.eye(2) * 0.1
C = np.array([[1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0]])

A_cl, B_cl, K, G = get_gains(params, x_eq, u_eq, Q, R, C)


# initial conditions
x0 = np.random.uniform(-3, 3, 6)

t_span = (0, 30)
t_eval = np.linspace(t_span[0], t_span[1], 600)

sol = solve_ivp(LQR_nonlinear_dynamics, t_span, x0, args=(K, G, params, x_eq, u_eq), t_eval=t_eval)

sol_linear = solve_ivp(LQR_linear_dynamics, t_span, x0, args=(K, G, A, B, x_eq, u_eq), t_eval=t_eval)

y = sol.y
y_linear = sol_linear.y

fig, ax = plt.subplots(6, 1, figsize=(8, 15))
state_names = [r'$x$', r'$y$', r'$\theta$', r'$v_x$', r'$v_y$', r'$\omega$']
color_names = niceplots.get_colors()
colors = [color_names['Blue'], color_names['Red'], color_names['Purple'], 
          color_names['Green'], color_names['Orange'], color_names['Navy']]
for i in range(6):
    niceplots.utils.adjust_spines(ax[i])

    ax[i].plot(sol.t, y[i, :], color=colors[i])
    ax[i].plot(sol_linear.t, y_linear[i, :], color=colors[i], linestyle='--')
    ax[i].set_ylabel(state_names[i], rotation=0)

    ax[i].yaxis.set_label_coords(-0.1, 0.5)

ax[-1].set_xlabel('Time (s)')

fig, ax = plt.subplots()
niceplots.utils.adjust_spines(ax)
ax.plot(y[0, :], y[1, :], label='Quadrotor Trajectory')
ax.plot(y_linear[0, :], y_linear[1, :], label='Linearized Trajectory', linestyle='--')
r = np.array([get_ref(ti) for ti in t_eval]).T
ax.plot(r[0, :], r[1, :], label='Reference Trajectory', linestyle='--')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$', rotation=0)
ax.legend(frameon=False)


plt.show()