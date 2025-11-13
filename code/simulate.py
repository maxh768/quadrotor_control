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

# optimal control weighting
Q = np.eye(6)
R = np.eye(2) * 0.1
C = np.array([[1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0]])
A_cl, B_cl, K, G = get_gains(params, x_eq, u_eq, Q, R, np.eye(6))
B_cl = 0

# initial conditions
x0 = np.random.uniform(-0.5, 0.5, 6)

t_span = (0, 30)
t_eval = np.linspace(t_span[0], t_span[1], 300)

sol = solve_ivp(LQR_nonlinear_dynamics, t_span, x0, args=(K, params, x_eq, u_eq), t_eval=t_eval)

y = sol.y

fig, ax = plt.subplots(6, 1, figsize=(8, 15))
state_names = [r'$x$', r'$y$', r'$\theta$', r'$v_x$', r'$v_y$', r'$\omega$']
color_names = niceplots.get_colors()
colors = [color_names['Blue'], color_names['Red'], color_names['Purple'], 
          color_names['Green'], color_names['Orange'], color_names['Navy']]
for i in range(6):
    niceplots.utils.adjust_spines(ax[i])

    ax[i].plot(sol.t, y[i, :], color=colors[i])
    ax[i].set_ylabel(state_names[i], rotation=0)

    ax[i].yaxis.set_label_coords(-0.1, 0.5)

ax[-1].set_xlabel('Time (s)')
plt.show()