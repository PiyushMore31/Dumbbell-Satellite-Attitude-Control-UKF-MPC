# Dumbbell Satellite Attitude Control — UKF + MPC

## Overview
Closed-loop attitude estimation and control system for a rigid dumbbell-shaped satellite on an eccentric orbit (e=0.6). An Unscented Kalman Filter (UKF) fuses gyroscope, sun sensor, and magnetometer measurements for full-state attitude estimation, which feeds a Model Predictive Control (MPC) algorithm that computes optimal reaction wheel torques over a 10-step prediction horizon subject to actuator constraints.

Developed as part of AAE 568 — Applied Optimal Control and Estimation, Purdue University (Spring 2026).

## System Description
- **Satellite:** Dumbbell configuration — two spherical masses connected by a thin rod
- **Orbit:** Eccentric (e=0.6), inclination 70°, SMA 7500 km
- **Attitude representation:** Unit quaternion
- **Sensors:** Gyroscope (0.01°/s noise) + Sun sensor (0.5° noise) + Magnetometer (1.0° noise) with first-order Gauss-Markov bias model

## Key Components

### UKF — Unscented Kalman Filter
- 6-state error-state formulation: attitude error (3) + gyroscope bias (3)
- Sigma point propagation using quaternion kinematics
- Measurement update fusing sun vector and magnetometer observations
- Gyroscope bias estimation with first-order Markov correlation model (τ = 2000s)

### MPC — Model Predictive Control
- Prediction horizon: N = 10 steps, sampling period Ts = 1s
- Linearized discrete-time error-state model about sun-pointing equilibrium
- Actuator torque constraints: |τ| ≤ 5×10⁻⁴ N·m per axis
- Quadratic cost: Q = 0.05I₆, R = 50I₃
- Condensed batch QP solved via MATLAB mpcActiveSetSolver

### Closed-Loop Architecture
UKF estimate → MPC → Plant (RK4) → UKF (receding horizon)

## Results
- MPC drives 65° initial misalignment to within 1° of sun-pointing in under 0.3 hours
- All actuator torques remain within ±5×10⁻⁴ N·m saturation limits
- UKF attitude errors settle within ±0.5° and remain bounded by 3σ envelopes for the full orbital period

## Dependencies
- MATLAB R2021a or later
- Aerospace Toolbox
- Model Predictive Control Toolbox

## Files
- `main.m` — Main simulation loop (UKF + MPC integration)
- `DARE_KF_Report.pdf` — Full technical report

## Course
AAE 568 — Applied Optimal Control and Estimation
Purdue University, Spring 2026
