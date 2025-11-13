import numpy as np
from dynamics import *
from controller import get_gains
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import niceplots

plt.style.use(niceplots.get_style('james-light'))


def add_direction_arrows(ax, x, y, num_arrows=6, color='k'):
    """Add arrows along a trajectory to show direction of travel."""
    if len(x) < 2 or num_arrows <= 0:
        return
    n = len(x)
    indices = np.linspace(0, n - 2, num_arrows, dtype=int)
    for idx in np.unique(indices):
        ax.annotate('', xy=(x[idx + 1], y[idx + 1]), xytext=(x[idx], y[idx]),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5))


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

A_cl, B_cl, K, G, L = get_gains(params, x_eq, u_eq, Q, R, C)


# initial conditions
x0 = np.random.uniform(-3, 3, 12)

t_span = (0, 15)
t_eval = np.linspace(t_span[0], t_span[1], 600)

sol = solve_ivp(closed_loop_system, t_span, x0, args=(params, A, B, C, K, G, L, x_eq, u_eq), t_eval=t_eval)

sol_linear = solve_ivp(closed_loop_system_linear, t_span, x0, args=(params, A, B, C, K, G, L, x_eq, u_eq), t_eval=t_eval)

y = sol.y[:6]

y_linear = sol_linear.y[:6]

fig, ax = plt.subplots(6, 1, figsize=(5, 8))
state_names = [r'$x$', r'$y$', r'$\theta$', r'$v_x$', r'$v_y$', r'$\omega$']
color_names = niceplots.get_colors()
colors = [color_names['Blue'], color_names['Red'], color_names['Purple'], 
          color_names['Green'], color_names['Orange'], color_names['Navy']]
for i in range(6):
    niceplots.utils.adjust_spines(ax[i])

    ax[i].plot(sol.t, y[i, :], color=colors[i], label='Nonlinear')
    ax[i].plot(sol_linear.t, y_linear[i, :], color=colors[i], linestyle='--', label='Linear')
    ax[i].set_ylabel(state_names[i], rotation=0)

    ax[i].yaxis.set_label_coords(-0.2, 0.5)

ax[-1].legend(frameon=False, loc='best')
ax[-1].set_xlabel('Time (s)')

niceplots.save_figs(fig, 'figures/state_trajectories.pdf', formats='pdf', bbox_inches='tight')

fig, ax = plt.subplots()
niceplots.utils.adjust_spines(ax)
traj_colors = {
    'nonlinear': color_names['Blue'],
    'linear': color_names['Red'],
    'ref': color_names['Green']
}

ax.plot(y[0, :], y[1, :], color=traj_colors['nonlinear'], label='Quadrotor Trajectory')
add_direction_arrows(ax, y[0, :], y[1, :], num_arrows=12, color=traj_colors['nonlinear'])
ax.plot(y_linear[0, :], y_linear[1, :], color=traj_colors['linear'], linestyle='--', label='Linearized Trajectory')
add_direction_arrows(ax, y_linear[0, :], y_linear[1, :], num_arrows=12, color=traj_colors['linear'])
r = np.array([get_ref(ti) for ti in t_eval]).T
ax.plot(r[0, :], r[1, :], color=traj_colors['ref'], linestyle=':', label='Reference Trajectory')
add_direction_arrows(ax, r[0, :], r[1, :], num_arrows=3, color=traj_colors['ref'])
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$', rotation=0)
ax.legend(frameon=False)


niceplots.save_figs(fig, 'figures/2D_trajectory.pdf', formats='pdf', bbox_inches='tight')

