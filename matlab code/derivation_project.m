%%
clear; clc;
%%
set(0, 'DefaultAxesFontName', 'Arial');
set(0,'DefaultAxesTitleFontWeight','normal');
set(groot,'defaultAxesFontSize',16);
%%

syms Tl Tr theta m g I l beta x y xdot ydot thetadot xddot yddot thetaddot

% define velocity norm
v_norm = sqrt(xdot^2 + ydot^2);

gamma = acos(xdot / v_norm);

force_x = -(Tl + Tr)*sin(theta);
force_y = (Tl + Tr)*cos(theta);

drag_mag = beta * v_norm^2;
drag_x = -drag_mag * cos(gamma);
drag_y = -drag_mag * sin(gamma);

tor_left = -Tl * l;
tor_right = Tr *l;

% x equation:
x_eq = m * xddot == force_x + drag_x;
simplify(x_eq);

% y equation:
y_eq = m * yddot == force_y + drag_y - m*g;
simplify(y_eq);

% theta eq:
theta_eq = I * thetaddot == tor_left + tor_right;
simplify(theta_eq);

%% linearization

% define states and inputs
q = [x;y;theta];

dq = [xdot;ydot;thetadot];

ddq = [xddot;yddot;thetaddot];

u = [Tl;Tr];

% get f = 0

f_x = lhs(x_eq) - rhs(x_eq);
f_y = lhs(y_eq) - rhs(y_eq);
f_theta = lhs(theta_eq) - rhs(theta_eq);

f = [f_x;f_y;f_theta];

% get dynamics in matrix form
% interia matrix: M = pf / p ddq
M = [diff(f,ddq(1)), diff(f,ddq(2)), diff(f,ddq(3))];
% C : pf / pdq
C = [diff(f,dq(1)), diff(f, dq(2)), diff(f, dq(3))];
% K : pf/pq
K = [diff(f, q(1)), diff(f, q(2)), diff(f, q(3))];
% N : pf/pu
N = [diff(f,u(1)), diff(f,u(2))];

% operating point
x_op = 1e-6;
y_op = 1e-6;
theta_op = 1e-6;
dx_op = 1e-6;
dy_op = 1e-6;
dtheta_op = 0.001;

Tl_op = m*g / 2;
Tr_op = m*g / 2;

% eval
M = subs(M,{x,y,theta,xdot,ydot,thetadot, Tl, Tr},{x_op, y_op, theta_op, dx_op, dy_op, dtheta_op, Tl_op, Tr_op});
C = subs(C,{x,y,theta,xdot,ydot,thetadot, Tl, Tr},{x_op, y_op, theta_op, dx_op, dy_op, dtheta_op, Tl_op, Tr_op});
K = subs(K,{x,y,theta,xdot,ydot,thetadot, Tl, Tr},{x_op, y_op, theta_op, dx_op, dy_op, dtheta_op, Tl_op, Tr_op});
N = subs(N,{x,y,theta,xdot,ydot,thetadot, Tl, Tr},{x_op, y_op, theta_op, dx_op, dy_op, dtheta_op, Tl_op, Tr_op});

% state space representation
A = [zeros(3,3), eye(3) ; -inv(M)*K, -inv(M)*C];
B = [zeros(3,2); -inv(M)*N];
A = simplify(A);
B = simplify(B);

% now define constants
m_use = 1.4;
I_use = 0.0211;
l_use = 0.159;
beta_use = 0.1365;
g_use = 9.81;

A = double(subs(A, {m, I, l, beta, g}, {m_use, I_use, l_use, beta_use, g_use}));
B = double(subs(B, {m, I, l, beta, g}, {m_use, I_use, l_use, beta_use, g_use}));

% stability of closed loop system
eig = eig(A);
eig;

% controllability:
P = ctrb(A, B);
rank(P);

% observability
% first define ouput y = [x, y] (don't include theta for now)
C = [1 0 0 0 0 0; 0 1 0 0 0 0]; %; 0 0 1 0 0 0];

O = obsv(A, C);
rank(O);

% transfer function
D = zeros(2);
sys_ss = ss(A, B, C, D);
sys_tf = tf(sys_ss)

% design LQR controller
Q = eye(6);
R = eye(2) * 0.2;

[K, s, p_c] = lqr(A, B, Q, R);

% Kalman Filter Design
Q = eye(6) * 0.01;
R = eye(2) * 0.01;
[L_t, fill, p_o] = lqr(transpose(A), transpose(C), Q, R);
L = transpose(L_t);

figure;
plot(real(p_o), imag(p_o), 'ro', 'MarkerSize', 12, 'DisplayName', 'Observer Poles'); hold on;
plot(real(p_c), imag(p_c), 'bx', 'MarkerSize', 12, 'DisplayName', 'Controller Poles');
xlabel('Re(p)');
ylabel('Im(p)', Rotation=0);
title('Observer and Controller Poles');
grid on;
legend('show', 'box', 'off', 'Location','south');
hold off;


% tracker with derivative tracking
G1 = inv(C * inv(B*K - A) * B);

G2 = ((C * ((B*K - A)^-2) * B)) / ( (C * ((B*K - A)^-1)*B)^2 );

sens_power = [1e-6;1e-6];
