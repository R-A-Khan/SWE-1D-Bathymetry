function [t,y] = RK34_FW_bath(dt, u0, b, tRange, trend, h)
% Usage: [t,y] = RK34_FW(dt, u0, tRange, trend, h)
%
% Explicit strongly stability preserving four stage third order Runge-Kutta
% Spiteri and Ruuth (SIAM J. Numer. Anal., 40(2): 469-491, 2002).
% Should be stable for CFL<=2 for a hyperbolic conservation law.
% 
% Input: 
% dt = time step size
% u0 = initial condition vector
% b  = bottom bathmetry 
% trend = approximation function
% 
% Output
% t = time values
% y = solution

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
numSteps = abs((tRange(2) - tRange(1))/dt);

t(1) = tRange(1);

for k = 1 : numSteps
    t(1,k+1) = t(1,k) + dt;
    
    U0 = y(:,k);
    U1 =  alpha(1,1)*U0                 + dt * beta(1,1)*trend(t,U0,h,b);
    U2 =  alpha(2,2)*U1                 + dt * beta(2,2)*trend(t,U1,h,b);
    U3 =  alpha(3,1)*U0 + alpha(3,3)*U2 + dt * beta(3,3)*trend(t,U2,h,b);
    U4 =  alpha(4,4)*U3                 + dt * beta(4,4)*trend(t,U3,h,b);
    y(:,k+1) = U4;
end