function dy = cent_diff_h(h, y)
% Usage: dy = cent_diff_h(h, y)
%
% Calculates central difference approximation to first derivative
% at mid-points with result at height points. 
% Accounts for periodic boundary conditions.
% Input: 
% y = data
% h = grid spacing
% Output: 
% dy = forward difference derivative at intermediate mid-points

N = length(y);
dy = zeros(N,1);

% Periodic boundary conditions
dy(N) = (y(1) - y(N))/h;
dy(1) = (y(2) - y(1))/h;

for i = 2:N-1
    dy(i) = (y(i+1) - y(i))/h;
end