function [eta_all, u_all, a, b, T, Y] = FW_solve_bath(h, dt, u0, b, tmin, tmax, x0_inds, trend)
% Usage:
% [eta_all, u_all, a, b, T, Y] = FW_solve(h, dt, u0, b, tmin, tmax, x0_inds, trend)
%
% Integrates SWE equations specified by function handle "trend" from tmin
% to tmax with time step dt using RK34.
%
% Input:
% N       = number of spatial grid points
% h       = grid spacing
% dt      = time step
% u0      = initial conditions for height and velocity
% b       = shape of bottom bathymetry
% tmin    = minimum time
% tmax    = maximum time (control time)
% x0_inds = indices of observation points
% trend   = function handle for trend (rhs of ode in time)
% 
% Output:
% eta_all = height at all positions and times
% u_all   = velocity at all positions and times
% a = height at each observation point x0 for all times
% b = height at control time tmax for all x
% T = numSteps x 1 vector of discrete time steps between tmin and tmax
% Y = 2N x numSteps matrix
%     1:N rows = height, N+1 : 2N rows = velocity
%     Rows   = values of height/velocity at point x_i for all times
%     Colums = values of height/velocity at time t_i for all positions

tRange = [tmin tmax];
[T,Y] = RK34_FW_bath(dt, u0, b, tRange, trend, h);

a = zeros(length(x0_inds),length(T));
for i= 1:length(x0_inds)
    a(i,:) = Y(x0_inds(i),:);
end
N = length(u0)/2;
b = Y(1:N,length(T));

% Height and velocity at all t and x
eta_all = Y(1:N, :);
u_all   = Y(N+1:2*N, :);
