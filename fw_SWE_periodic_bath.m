function trend = fw_SWE_periodic_bath(t, y, h, b)
% Usage: trend = fw_SWE_periodic(t, y, h)
% Finds trend of linear SWE equations with periodic boundary conditions
% using second order finite difference / finite volume approximation
% grid size is 2*N, with dh/dt evaluated at centres and du/dt at edges
% use mid2edge on bathymetry(at centres) so u*(1-b) is all on edge points
%
% 1 to N rows of trend represent dh/dt = -d(u(1-b))/dx 
% = -1/h * ([u*(1-b)]_i+1 - [u*(1-b)]_i)
% N+1 to 2*N rows of trend represent du/dt = - dh/dx = -1/h * (h_i - h_i-1)
%
% Input:
% t = time (not used)
% y = input variables (height and velocity in one vector)
% h = spacing of grid points
% Output:
% trend = trend for both height and velocity equations in one vector

N = length(y)/2;
eta = y(1:N);
u   = y(N+1:2*N);

trend = zeros(2*N, 1);
trend(1:N)     =  -cent_diff_h(h, mid2edge(1+eta-b).*u);
trend(N+1:2*N) =  -cent_diff_u(h, edge2mid(u.^2/2) + eta);
