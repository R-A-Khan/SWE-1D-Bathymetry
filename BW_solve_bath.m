function [a, b, T, Y] = BW_solve_bath(h, dt, u0, tmin, tmax, x0_inds, trend, obs, eta, beta, u_all, eta_all, T_obs, T_eta)
% Usage: 
% [b, T, Y] = BW_solve(h, dt, u0, tmin, tmax, x0_inds, trend, obs, eta, u_all, eta_all, T_obs, T_eta, nl)
%
% Solves backwards (adjoint) equations using observations results from
% forward solver
%
% Input:
% h       = grid spacing
% dt      = time step
% u0      = initial conditions for height and velocity
% tmin    = minimum time
% tmax    = maximum time (control time)
% x0_inds = indices of observation positions
% trend   = function handle for trend (rhs of ode in time)
% obs     = vector of height at observation position from t=tmin to t=tmax
% eta     = height at control time for all positions
% beta    = bathymetry at control time all positions
% u_all   = velocity at all positions and times
% eta_all = height at all positions and times
% T_obs   = times for all observations
% T_eta   = times for all height data from forward solver
% nl      = true = nonlinear equations, false = linear equations
%
% Output:
% a = velocity at time tmin for all positions
% b = height at time tmin for all positions
% T = numSteps x 1 vector of discrete time steps between tmin and tmax
% Y = 2N x numSteps matrix
%     1:N rows = height, N+1 : 2N rows = velocity
%     Rows   = values of height/velocity at point x_i for all times
%     Colums = values of height/velocity at time t_i for all positions

% Solving backwards starting t=tmax to t=tmin
% Flip eta and obs e.g. obs_1 corresponds to t=tmax


obs = fliplr(obs);
eta = fliplr(eta);
beta = fliplr(beta);
u_all   = fliplr(u_all);
eta_all = fliplr(eta_all);

tRange = [tmin tmax];
[T,Y] = RK34_BW_bath(dt, u0, tRange, trend, h, x0_inds, obs, eta, beta, u_all, eta_all);

N = length(u0)/2;
a = Y(N+1:2*N,length(T));
b = Y(1:N,length(T));

