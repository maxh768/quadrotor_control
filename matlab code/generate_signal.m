
t = 0:0.01:120;
radius = 0.2;
period = 30.0;
omega = 2 * pi / period;
x_ref = radius * cos(omega * t);
y_ref = radius * sin(omega * t);

dx_ref = - radius * omega * sin(omega*t);
dy_ref = radius*omega*cos(omega*t);

reference_traj = [x_ref', y_ref'];
reference_traj = timeseries(reference_traj, t);

ref_traf_der = [dx_ref', dy_ref'];
ref_traj_der = timeseries(ref_traf_der, t);