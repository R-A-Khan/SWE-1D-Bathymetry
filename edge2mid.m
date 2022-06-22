function u_centre = edge2mid(u_edge)
% Usage: u_centre = edge2mid(u_edge)
% Linear interpolation (second order) from cell edge (velocity) points to cell centre (height)
% points. Assumes doubly periodic boundary conditions and uniform grid
% spacing.
%
% Input:
% u_edge = values at cell edges
% Output:
% u_centre = values at cell centres

N = length(u_edge);
u_centre = zeros(N,1);

if N >1 

u_centre(1) = (u_edge(1) + u_edge(2))/2;
u_centre(N) = (u_edge(N) + u_edge(1))/2;

for i = 2:N-1
    u_centre(i) = (u_edge(i) + u_edge(i+1))/2;
end
else u_centre = u_edge;
end