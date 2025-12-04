%%
set(0, 'DefaultAxesFontName', 'Arial');
set(0,'DefaultAxesTitleFontWeight','normal');
set(groot,'defaultAxesFontSize',16);
%%

t = full_state_true.time;

% nonlinear states
states_true = full_state_true.Data;
states_linear = full_state_linear.Data;

x_true = states_true(1,:);
y_true = states_true(2,:);
theta_true = states_true(3,:);
dx_true = states_true(4,:);
dy_true = states_true(5,:);
dtheta_true = states_true(6,:);

x_linear = states_linear(:,1);
y_linear = states_linear(:,2);
theta_linear = states_linear(:,3);
dx_linear = states_linear(:,4);
dy_linear = states_linear(:,5);
dtheta_linear = states_linear(:,6);

u_hat = m_use * g_use / 2;

% nonlinear control input
u_true = ctrl_input_true.Data;
u_linear = ctrl_input_linear.Data;
u_linear = u_linear + u_hat;
Tl_true = u_true(1,:);
Tr_true = u_true(2,:);
Tl_linear = u_linear(:,1);
Tr_linear = u_linear(:,2);



% nonlinear observer state
obs_true = obs_state_true.Data;
obs_linear = obs_state_linear.Data;
obs_x_true = obs_true(1,:);
obs_y_true = obs_true(2,:);
obs_theta_true = obs_true(3,:);
obs_dx_true = obs_true(4,:);
obs_dy_true = obs_true(5,:);
obs_dtheta_true = obs_true(6,:);

obs_x_linear = obs_linear(:,1);
obs_y_linear = obs_linear(:,2);
obs_theta_linear = obs_linear(:,3);
obs_dx_linear = obs_linear(:,4);
obs_dy_linear = obs_linear(:,5);
obs_dtheta_linear = obs_linear(:,6);

% referemce signal
ref = ref_sig_true.Data;
ref_t = ref_sig_true.time;
x_ref = ref(:,1);
y_ref = ref(:,2);

% state plots
figure;

% True states
subplot(6,2,1)
plot(t, x_true, 'Color', 'b', 'LineWidth', 2)
title('Nonlinear System (True) States')
ylabel('x', 'Rotation', 0)
grid on

subplot(6,2,3)
plot(t, y_true, 'Color', 'b', 'LineWidth', 2)
ylabel('y', 'Rotation', 0)
grid on

subplot(6,2,5)
plot(t, theta_true, 'Color', 'b', 'LineWidth', 2)
ylabel('\theta', 'Rotation', 0)
grid on

subplot(6,2,7)
plot(t, dx_true, 'Color', 'b', 'LineWidth', 2)
ylabel('dx', 'Rotation', 0)
grid on

subplot(6,2,9)
plot(t, dy_true, 'Color', 'b', 'LineWidth', 2)
ylabel('dy', 'Rotation', 0)
grid on

subplot(6,2,11)
plot(t, dtheta_true, 'Color', 'b', 'LineWidth', 2)
ylabel('d\theta', 'Rotation', 0)
grid on
xlabel('Time (s)')

% Linear states
subplot(6,2,2)
plot(t, x_linear, 'Color', 'r', 'LineWidth', 2)
title('Linear System States')
grid on

subplot(6,2,4)
plot(t, y_linear, 'Color', 'r', 'LineWidth', 2)
grid on

subplot(6,2,6)
plot(t, theta_linear, 'Color', 'r', 'LineWidth', 2)
grid on

subplot(6,2,8)
plot(t, dx_linear, 'Color', 'r', 'LineWidth', 2)
grid on

subplot(6,2,10)
plot(t, dy_linear, 'Color', 'r', 'LineWidth', 2)
grid on

subplot(6,2,12)
plot(t, dtheta_linear, 'Color', 'r', 'LineWidth', 2)
grid on
xlabel('Time (s)')

% control plots
figure;

subplot(1,2,1)
plot(t, Tl_true, 'LineWidth', 2)
hold on
plot(t, Tr_true, 'LineWidth', 2)
grid on
xlabel('Time (s)')
title('Nonlinear System Control Inputs')
ylabel('u', 'Rotation', 0)
legend('Tl', 'Tr', 'Box', 'off', 'Location','northwest')
ylim([0, 10])

subplot(1,2,2)
plot(t, Tl_linear, 'LineWidth', 2)
hold on
plot(t, Tr_linear, 'LineWidth', 2)
grid on
ylim([-30 30])
title('Linear System Control Inputs')
xlabel('Time (s)')
ylim([0, 10])


% Error calculation for nonlinear system
error_x_true = obs_x_true - x_true;
error_y_true = obs_y_true - y_true;
error_theta_true = obs_theta_true - theta_true;
error_dx_true = obs_dx_true - dx_true;
error_dy_true = obs_dy_true - dy_true;
error_dtheta_true = obs_dtheta_true - dtheta_true;

norm_error_true = sqrt(error_x_true.^2 + error_y_true.^2 + error_theta_true.^2 + ...
                        error_dx_true.^2 + error_dy_true.^2 + error_dtheta_true.^2);

% Error calculation for linear system
error_x_linear = obs_x_linear - x_linear;
error_y_linear = obs_y_linear - y_linear;
error_theta_linear = obs_theta_linear - theta_linear;
error_dx_linear = obs_dx_linear - dx_linear;
error_dy_linear = obs_dy_linear - dy_linear;
error_dtheta_linear = obs_dtheta_linear - dtheta_linear;

norm_error_linear = sqrt(error_x_linear.^2 + error_y_linear.^2 + error_theta_linear.^2 + ...
                          error_dx_linear.^2 + error_dy_linear.^2 + error_dtheta_linear.^2);

% Plotting the norm errors
figure;

% Nonlinear system error
subplot(1,2,1)
plot(t, norm_error_true, 'Color', 'b', 'LineWidth', 2)
title('Nonlinear System Observer Error')
xlabel('Time (s)')
ylabel('Norm Error', Rotation=0)
grid on

% Linear system error
subplot(1,2,2)
plot(t, norm_error_linear, 'Color', 'r', 'LineWidth', 2)
ylim([0, 5])
title('Linear System Observer Error')
xlabel('Time (s)')
grid on

figure;
plot(x_ref, y_ref, 'Color', 'k', 'LineWidth', 2)
hold on
plot(x_linear, y_linear, 'Color', 'r', 'LineStyle', '--')
plot(x_true, y_true, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2)
title(['Observer + FSFB Trajectories with period T = ', num2str(period), ' s'])
xlabel('x')
ylabel('y', 'Rotation', 0)
legend('Reference', 'Linear', 'True', 'Location', 'northeast', 'Box', 'off')


