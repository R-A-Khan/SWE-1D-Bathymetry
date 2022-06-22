function h_edge = mid2edge(h_centre)
% Usage: h_edge = mid2edge(h_centre)
% Linear interpolation (second order) from cell centres (height) points to cell edges (velocity)
% points. Assumes doubly periodic boundary conditions and uniform grid
% spacing.
%
% Input:
% h_centre = values at cell centres
% Output:
% h_edge = values at cell edges

N = length(h_centre);
h_edge = zeros(N,1);

h_edge(1) = (h_centre(N) + h_centre(1))/2;
h_edge(N) = (h_centre(N-1) + h_centre(N))/2;

for i = 2:N-1
    h_edge(i) = (h_centre(i-1) + h_centre(i))/2;
end

