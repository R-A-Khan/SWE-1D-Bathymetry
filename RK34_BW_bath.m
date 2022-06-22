function [t,y] = RK34_BW_bath(dt, u0, tRange, trend, h, x0, obs, eta, bath, u_all, eta_all)
% Usage:
% [t,y] = RK34_BW(dt, u0, tRange, trend, h, x0, obs, eta, u_all, eta_all, T_obs, T_eta)
%
% Input: 
% dt = time step size
% u0 = initial condition vector
% tRange = range of time values
% trend = approximation function
%
% Output:
% t = time values
% y = solution

% Explicit strongly stability preserving four stage third order Runge-Kutta
% Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002).
% Should be stable for CFL<=2 for a hyperbolic conservation law.

% obs  = I x length(T)  vector for height at x0 points for all t given
% actual bathymetry
% eta  = I x length(T)  vector for height at x0 points for all t given
% "guess" bathymetry
% beta = I x length(T)  vector for bathymetry at x0 points for all t given
% "guess" bathymetry
% Tx1 vectors from t_0 to t_T


alpha = [ 1   0  0  0 ; ...
          0   1  0  0 ; ...
         2/3  0 1/3 0 ; ...
          0   0  0  1 ];

beta = [1/2  0   0    0 ; ...
         0  1/2  0    0 ; ...
         0   0  1/6   0 ; ...
         0   0   0   1/2 ];

% Setting the initial condition to the first column of y
y(:,1) = u0;

% Finding numSteps
numSteps = abs(( tRange(2) - tRange(1) ) / dt);

t(1,1) = tRange(1);
for k = 1 :  numSteps
    t(1,k+1) = t(1,k) + dt;
    
    U0 = y(:,k);
    U1 =  alpha(1,1)*U0                 + dt * beta(1,1)*trend(t,U0,h,x0,obs(:,k),eta(:,k), bath, u_all(:,k), eta_all(:,k));
    U2 =  alpha(2,2)*U1                 + dt * beta(2,2)*trend(t,U1,h,x0,obs(:,k),eta(:,k), bath, u_all(:,k), eta_all(:,k));
    U3 =  alpha(3,1)*U0 + alpha(3,3)*U2 + dt * beta(3,3)*trend(t,U2,h,x0,obs(:,k),eta(:,k), bath, u_all(:,k), eta_all(:,k));
    U4 =  alpha(4,4)*U3                 + dt * beta(4,4)*trend(t,U3,h,x0,obs(:,k),eta(:,k), bath, u_all(:,k), eta_all(:,k));
    y(:,k+1) = U4;
end

