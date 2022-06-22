function dy = cent_diff_u(h, y)
% Usage: dy = cent_diff_h(h, y)
%
% Calculates central difference approximation to first derivative
% at mid-points with result at velocity points. 
% Accounts for periodic boundary conditions.
% Input: 
% y = data
% h = grid spacing
% Output: 
% dy = backward difference derivative at intermediate mid-points

N = length(y);
dy = zeros(N,1);

% Periodic boundary conditions
dy(N) = (y(N) - y(N-1))/h;
dy(1) = (y(1) - y(N))/h;

for i = 2:N-1
    dy(i) = (y(i) - y(i-1))/h;
end