% AAE 568 - Attitude Control: UKF + MPC
% Project DARE - Attitude Control: UKF + MPC

clc; clear all; close all;
format long g

%% Physical Parameters
mu  = 398600 * 1000^3;     % Earth grav param [m^3/s^2]
SMA = 7500 * 1000;         % semi-major axis [m]

% Dumbbell satellite geometry
ms  = 5;   % sphere mass [kg]
Rsp = 0.1;   % sphere radius [m]
l   = 1.0;   % center-to-center distance [m]
mr  = 1;     % rod mass [kg]

I = dumbbellInertia(ms, Rsp, l, mr);   % 3x3 inertia tensor

%% Orbit Setup
e    = 0.6;    % eccentricity
TA0  = 45;     % initial true anomaly [deg]
w    = 10;     % argument of periapsis [deg]
inc  = 70;     % inclination [deg]
RAAN = 70;     % right ascension of ascending node [deg]

T     = 2*pi*sqrt(SMA^3/mu);   % orbit period [s]
tf    = 1*T;
dt    = 1.0;
tspan = 0:dt:tf;
N     = length(tspan);

%% Initial Conditions
% Body axis Euler angles [deg]
Roll = 80;  Pitch = 120;  Yaw = 160;

% Initial angular rates [deg/s]
wx = 0.5;  wy = 0.3;  wz = 0.2;

% Convert to quaternion (scalar-first [q0 q1 q2 q3])
eul0      = deg2rad([Yaw, Pitch, Roll]);
q_initial = eul2quat(eul0, 'ZYX');

% Initial position/velocity from orbital elements
r0      = SMA*(1-e^2)/(1+e*cosd(TA0));
h       = sqrt(mu*SMA*(1-e^2));

r0_peri = r0*[cosd(TA0); sind(TA0); 0];
v0_peri = (mu/h*e*sind(TA0))*[cosd(TA0); sind(TA0); 0] + ...
          (h/r0)*[-sind(TA0); cosd(TA0); 0];

% Perifocal -> ECI rotation
R3_W = [cosd(RAAN) -sind(RAAN) 0; sind(RAAN) cosd(RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cosd(inc) -sind(inc); 0 sind(inc) cosd(inc)];
R3_w = [cosd(w) -sind(w) 0; sind(w) cosd(w) 0; 0 0 1];

r0_eci = R3_W*R1_i*R3_w*r0_peri;
v0_eci = R3_W*R1_i*R3_w*v0_peri;

omega0 = deg2rad([wx; wy; wz]);
state0 = [r0_eci; v0_eci; q_initial'; omega0];   % 13x1
state  = state0;
f_dyn  = @(t,x,tau) augEoM(t, x, mu, I, tau);

%% Sensor Parameters
% Reference vectors in inertial frame
sun_I = [1; 0; 0];    sun_I = sun_I/norm(sun_I);
mag_I = [0.2; 0.3; -0.9];  mag_I = mag_I/norm(mag_I);

% Measurement noise
sigma_sun  = deg2rad(0.5);
sigma_mag  = deg2rad(1.0);
sigma_gyro = deg2rad(0.01);

% Gyro bias model (first-order Markov)
bias_stab_std    = deg2rad(0.5)/3600;
correlation_time = 2000;
alpha_bias       = exp(-dt/correlation_time);
sig_torque       = 1e-6;

%% UKF Initialization
n_s      = 6;          % error-state: [dtheta(3); dbias(3)]
alpha_u  = 1;
beta_u   = 2;
kappa_u  = 0;
lambda_u = alpha_u^2*(n_s + kappa_u) - n_s;
gamma_u  = sqrt(n_s + lambda_u);

% Sigma point weights
Wm    = [lambda_u/(n_s+lambda_u), repmat(1/(2*(n_s+lambda_u)), 1, 2*n_s)];
Wc    = Wm;
Wc(1) = Wc(1) + (1 - alpha_u^2 + beta_u);

% Initial uncertainties
sigma_init      = deg2rad(5);
sigma_bias_init = deg2rad(30)/3600;
sigma_bias_proc = sqrt(1-alpha_bias^2)*bias_stab_std;

% Noise covs for UKF
Q_ukf  = blkdiag((sigma_gyro*dt*1.3)^2*eye(3), (sigma_bias_proc^2*1.5)*eye(3));
R_meas = blkdiag(sigma_sun^2*eye(3), sigma_mag^2*eye(3));
P_ukf  = blkdiag((sigma_init^2*7)*eye(3), (sigma_bias_init^2)*eye(3));

% perturb init est
delta_theta0 = sigma_init*randn(3,1);
dq0          = smallAngleQuat(delta_theta0);
qhat         = quatmultiply(q_initial, dq0')';
qhat         = qhat/norm(qhat);

bias_hat = zeros(3,1);
bias     = deg2rad(2)/3600*randn(3,1);
gyro_k   = omega0 + bias + sigma_gyro*randn(3,1);

%% MPC Setup
Ts    = 1.0;
Ad    = [eye(3), (Ts/2)*eye(3); zeros(3,3), eye(3)];
Bd    = [zeros(3,3); Ts*inv(I)];
A_mpc = Ad;
B_mpc = Bd;

n_mpc   = 6;
m_mpc   = 3;
N_hor   = 10;
tau_max = 0.0005;

Q_mpc  = 0.05*eye(6);
R_mpc  = 50*eye(3);
P_term = Q_mpc;
x_ref  = zeros(6,1);

% Build condensed QP matrices
Z_blk = zeros(n_mpc, m_mpc);

Qbar = blkdiag(Q_mpc,Q_mpc,Q_mpc,Q_mpc,Q_mpc,Q_mpc,Q_mpc,Q_mpc,Q_mpc,P_term);
Rbar = blkdiag(R_mpc,R_mpc,R_mpc,R_mpc,R_mpc,R_mpc,R_mpc,R_mpc,R_mpc,R_mpc);

% Prediction horizon matrices
H_mpc = [A_mpc; A_mpc^2; A_mpc^3; A_mpc^4; A_mpc^5;
         A_mpc^6; A_mpc^7; A_mpc^8; A_mpc^9; A_mpc^10];

G_mpc = [B_mpc        Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk;
         A_mpc*B_mpc  B_mpc    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk;
         A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk;
         A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk;
         A_mpc^4*B_mpc  A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk    Z_blk    Z_blk    Z_blk    Z_blk;
         A_mpc^5*B_mpc  A_mpc^4*B_mpc  A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk    Z_blk    Z_blk    Z_blk;
         A_mpc^6*B_mpc  A_mpc^5*B_mpc  A_mpc^4*B_mpc  A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk    Z_blk    Z_blk;
         A_mpc^7*B_mpc  A_mpc^6*B_mpc  A_mpc^5*B_mpc  A_mpc^4*B_mpc  A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk    Z_blk;
         A_mpc^8*B_mpc  A_mpc^7*B_mpc  A_mpc^6*B_mpc  A_mpc^5*B_mpc  A_mpc^4*B_mpc  A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc    Z_blk;
         A_mpc^9*B_mpc  A_mpc^8*B_mpc  A_mpc^7*B_mpc  A_mpc^6*B_mpc  A_mpc^5*B_mpc  A_mpc^4*B_mpc  A_mpc^3*B_mpc  A_mpc^2*B_mpc  A_mpc*B_mpc  B_mpc];


F_cost = G_mpc'*Qbar*H_mpc;
L_cost = G_mpc'*Qbar*G_mpc + Rbar;

E_con = [eye(m_mpc*N_hor); -eye(m_mpc*N_hor)];
W_con = tau_max*ones(2*m_mpc*N_hor, 1);

opt_mpc = mpcActiveSetOptions;
opt_mpc.IntegrityChecks = false;
iA = false(size(W_con));

%% Storage Allocation
state_hist    = zeros(13, N);
qhat_hist     = zeros(N, 4);
bias_hat_hist = zeros(N, 3);
P_hist        = zeros(6, 6, N);
X_store       = zeros(6, N);
U_store       = zeros(3, N-1);

state_hist(:,1)    = state0;
qhat_hist(1,:)     = qhat';
bias_hat_hist(1,:) = bias_hat';
P_hist(:,:,1)      = P_ukf;
X_store(:,1)       = [qhat(2:4); gyro_k - bias_hat];

%% Main Loop
for k = 1:N-1

    % --- MPC: compute control torque from current UKF estimate ---
    x_est = [qhat(2:4); gyro_k - bias_hat];
    f_vec = F_cost*(x_est - x_ref);
    [U_opt, ~, iA] = mpcActiveSetSolver(L_cost, f_vec, E_con, W_con, [], zeros(0,1), iA, opt_mpc);
    U_store(:,k) = U_opt(1:3);
    
    % --- Propagate true plant (RK4) ---
    tau_ctrl = U_store(:,k) + sig_torque*randn(3,1);
    k1 = f_dyn(tspan(k),        state,          tau_ctrl);
    k2 = f_dyn(tspan(k)+dt/2,   state+dt/2*k1,  tau_ctrl);
    k3 = f_dyn(tspan(k)+dt/2,   state+dt/2*k2,  tau_ctrl);
    k4 = f_dyn(tspan(k)+dt,     state+dt*k3,    tau_ctrl);
    state        = state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    state(7:10)  = state(7:10)/norm(state(7:10));   % re-normalize quaternion
    state_hist(:,k+1) = state;

    % --- Simulate sensors ---
    omega_true_k = state(11:13);
    bias   = alpha_bias*bias + sqrt(1-alpha_bias^2)*bias_stab_std*randn(3,1);
    gyro_k = omega_true_k + bias + sigma_gyro*randn(3,1);

    q_true_k = state(7:10);
    RbI_true = quat2dcm(q_true_k');
    sun_meas = quat2dcm(smallAngleQuat(sigma_sun*randn(3,1))') * RbI_true * sun_I;
    mag_meas = quat2dcm(smallAngleQuat(sigma_mag*randn(3,1))') * RbI_true * mag_I;
    sun_meas = sun_meas/norm(sun_meas);
    mag_meas = mag_meas/norm(mag_meas);

    % --- UKF: generate sigma points ---
    S_chol = chol(P_ukf + 1e-10*eye(n_s), 'lower');
    Xi = zeros(n_s, 2*n_s+1);
    for i = 1:n_s
        Xi(:, i+1)      =  gamma_u*S_chol(:,i);
        Xi(:, i+1+n_s)  = -gamma_u*S_chol(:,i);
    end

    q_sigma    = zeros(4, 2*n_s+1);
    bias_sigma = zeros(3, 2*n_s+1);
    for i = 1:2*n_s+1
        dq_s           = smallAngleQuat(Xi(1:3,i));
        q_sigma(:,i)   = quatmultiply(dq_s', qhat')';
        bias_sigma(:,i)= bias_hat + Xi(4:6,i);
    end

    % --- UKF: propagate sigma points ---
    q_prop    = zeros(4, 2*n_s+1);
    bias_prop = zeros(3, 2*n_s+1);
    for i = 1:2*n_s+1
        w_c = gyro_k - bias_sigma(:,i);
        Om  = [0 -w_c'; w_c -skew(w_c)];
        q_new        = q_sigma(:,i) + 0.5*Om*q_sigma(:,i)*dt;
        q_prop(:,i)    = q_new/norm(q_new);
        bias_prop(:,i) = alpha_bias*bias_sigma(:,i);
    end

    % --- UKF: weighted mean ---
    bias_hat = zeros(3,1);
    for i = 1:2*n_s+1
        bias_hat = bias_hat + Wm(i)*bias_prop(:,i);
    end

    q_ref2     = q_prop(:,1);
    dtheta_sum = zeros(3,1);
    for i = 1:2*n_s+1
        q_err = quatmultiply(q_prop(:,i)', quatinv(q_ref2'))';
        if q_err(1) < 0, q_err = -q_err; end
        dtheta_sum = dtheta_sum + Wm(i)*2*q_err(2:4);
    end
    dq_mean = smallAngleQuat(dtheta_sum);
    qhat    = quatmultiply(dq_mean', q_ref2')';
    qhat    = qhat/norm(qhat);

    % --- UKF: predicted covariance ---
    P_ukf = Q_ukf;
    for i = 1:2*n_s+1
        q_err  = quatmultiply(q_prop(:,i)', quatinv(qhat'))';
        if q_err(1) < 0, q_err = -q_err; end
        dtheta = 2*q_err(2:4);
        dbias  = bias_prop(:,i) - bias_hat;
        dx     = [dtheta; dbias];
        P_ukf  = P_ukf + Wc(i)*(dx*dx');
    end

    % --- UKF: predicted measurements ---
    Z_sig = zeros(6, 2*n_s+1);
    for i = 1:2*n_s+1
        RbI_i      = quat2dcm(q_prop(:,i)');
        sp         = RbI_i*sun_I;  sp = sp/norm(sp);
        mp         = RbI_i*mag_I;  mp = mp/norm(mp);
        Z_sig(:,i) = [sp; mp];
    end
    z_hat = zeros(6,1);
    for i = 1:2*n_s+1
        z_hat = z_hat + Wm(i)*Z_sig(:,i);
    end

    % --- UKF: innovation covariance ---
    S_cov = R_meas;
    Pxz   = zeros(n_s, 6);
    for i = 1:2*n_s+1
        dz    = Z_sig(:,i) - z_hat;
        q_err = quatmultiply(q_prop(:,i)', quatinv(qhat'))';
        if q_err(1) < 0, q_err = -q_err; end
        dtheta = 2*q_err(2:4);
        dbias  = bias_prop(:,i) - bias_hat;
        dx     = [dtheta; dbias];
        S_cov  = S_cov + Wc(i)*(dz*dz');
        Pxz    = Pxz   + Wc(i)*(dx*dz');
    end

    % --- UKF: measurement update ---
    K        = Pxz/S_cov;
    dx_upd   = K*([sun_meas; mag_meas] - z_hat);
    qhat     = quatmultiply(smallAngleQuat(dx_upd(1:3))', qhat')';
    qhat     = qhat/norm(qhat);
    bias_hat = bias_hat + dx_upd(4:6);
    P_ukf    = P_ukf - K*S_cov*K';
    P_ukf    = 0.5*(P_ukf + P_ukf') + 1e-10*eye(n_s);   % enforce symmetry

    % Store
    qhat_hist(k+1,:)     = qhat';
    bias_hat_hist(k+1,:) = bias_hat';
    P_hist(:,:,k+1)      = P_ukf;
    X_store(:,k+1)       = [qhat(2:4); gyro_k - bias_hat];

end

fprintf('Simulation complete: %d steps\n', N-1);

%% Post-Processing: Angle to Sun
t_plot        = (0:N-1)*dt/3600;
angle_sun_hist = zeros(N,1);
for k = 1:N
    q_k  = state_hist(7:10,k);
    R_bI = quat2dcm(q_k');
    bx   = R_bI'*[1;0;0];
    angle_sun_hist(k) = acosd(dot(bx, sun_I));
end

%% UKF Attitude Error
att_err = zeros(N,3);
for k = 1:N
    qt   = state_hist(7:10,k)/norm(state_hist(7:10,k));
    qe   = qhat_hist(k,:)'/norm(qhat_hist(k,:)');
    qerr = quatmultiply(qt', quatinv(qe'))';
    if qerr(1) < 0, qerr = -qerr; end
    att_err(k,:) = 2*qerr(2:4)';
end

sig_att = zeros(N,3);
for k = 1:N
    P_att      = P_hist(1:3,1:3,k);
    sig_att(k,:) = 3*sqrt(diag(P_att))';
end

t_zoom = [0.75*T/3600, T/3600];   % last 25% of orbit

%% Plots
figure
plot(t_plot, state_hist(7,:), 'LineWidth', 1.5, 'DisplayName', 'q_0')
hold on
plot(t_plot, state_hist(8,:), 'LineWidth', 1.5, 'DisplayName', 'q_1')
plot(t_plot, state_hist(9,:), 'LineWidth', 1.5, 'DisplayName', 'q_2')
plot(t_plot, state_hist(10,:), 'LineWidth', 1.5, 'DisplayName', 'q_3')
legend; grid on
title('True Quaternion States')
xlabel('Time (hrs)')
xlim(t_zoom)

figure
plot(t_plot, state_hist(11,:), 'LineWidth', 1.5, 'DisplayName', '\omega_x')
hold on
plot(t_plot, state_hist(12,:), 'LineWidth', 1.5, 'DisplayName', '\omega_y')
plot(t_plot, state_hist(13,:), 'LineWidth', 1.5, 'DisplayName', '\omega_z')
legend; grid on
title('True Angular Rates')
xlabel('Time (hrs)')
xlim(t_zoom)

figure
plot(t_plot(1:end-1), U_store', 'LineWidth', 1.5)
yline( tau_max, 'r--', 'LineWidth', 1.5)
yline(-tau_max, 'r--', 'LineWidth', 1.5)
legend('\tau_x', '\tau_y', '\tau_z')
title('Control Torques')
xlabel('Time (hrs)'); grid on
% xlim(t_zoom)
xlim([0 1.8]) 

figure
plot(t_plot, rad2deg(att_err(:,1)), 'r', 'DisplayName', '\delta\theta_x')
hold on
plot(t_plot, rad2deg(att_err(:,2)), 'g', 'DisplayName', '\delta\theta_y')
plot(t_plot, rad2deg(att_err(:,3)), 'b', 'DisplayName', '\delta\theta_z')
grid on; xlabel('Time (hrs)'); ylabel('Error (deg)')
title('UKF Attitude Estimation Error'); legend

figure
hold on
plot(t_plot, rad2deg(att_err(:,1)), 'r', 'DisplayName', '\delta\theta_x')
plot(t_plot, rad2deg(att_err(:,2)), 'g', 'DisplayName', '\delta\theta_y')
plot(t_plot, rad2deg(att_err(:,3)), 'b', 'DisplayName', '\delta\theta_z')
plot(t_plot,  rad2deg(sig_att(:,1)), 'r--', 'DisplayName', '3\sigma_x')
plot(t_plot, -rad2deg(sig_att(:,1)), 'r--', 'HandleVisibility', 'off')
plot(t_plot,  rad2deg(sig_att(:,2)), 'g--', 'DisplayName', '3\sigma_y')
plot(t_plot, -rad2deg(sig_att(:,2)), 'g--', 'HandleVisibility', 'off')
plot(t_plot,  rad2deg(sig_att(:,3)), 'b--', 'DisplayName', '3\sigma_z')
plot(t_plot, -rad2deg(sig_att(:,3)), 'b--', 'HandleVisibility', 'off')
xlabel('Time (hrs)'); ylabel('Attitude Error (deg)')
title('UKF Attitude Estimation Error with 3\sigma Bounds')
legend; grid on
ylim([-4 4])

figure
plot(t_plot, angle_sun_hist, 'b', 'LineWidth', 1.5)
yline(0, 'r--', 'LineWidth', 1.5)
title('Angle to Sun Vector')
xlabel('Time (hrs)'); ylabel('Angle (deg)')
legend('Angle to Sun', 'Reference (0 deg)'); grid on

%% Functions
function I = dumbbellInertia(ms, R, L, mr)
    Ixx = 2*(2/5)*ms*R^2;
    Iyy = 2*((2/5)*ms*R^2 + ms*(L/2)^2) + (1/12)*mr*(L-2*R)^2;
    Izz = Iyy;
    I   = diag([Ixx, Iyy, Izz]);
end

function xdot = augEoM(~, x, mu, I, tau)
    r = x(1:3);  v = x(4:6);
    q = x(7:10); q = q/norm(q);
    omega = x(11:13);
    Om    = [0 -omega'; omega -skew(omega)];
    xdot  = [v; -mu*r/norm(r)^3; 0.5*Om*q; I\(tau - cross(omega, I*omega))];
end

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function dq = smallAngleQuat(dtheta)
    dtheta = dtheta(:);
    dq = [1; 0.5*dtheta];
    dq = dq/norm(dq);
end

function qinv = quatinv(q)
    q = q(:)';
    qinv = [q(1), -q(2:4)];
end