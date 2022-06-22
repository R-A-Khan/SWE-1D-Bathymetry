function  [x0_inds, x0_pts] = mult_x0(X,x0_vals)
%
% Usage: 
% [x0_inds, x0_pts] = mult_x0(X,x0_vals)
%
% Given desired positions of observation points and position coordinates of 
% grid points returns grid indices of observation points and height at
% observation points
%
% Input:
% X = position coordinates of grid points
% x0_vals = positions of observation points (real number)
%
% Output:
% x0_inds = indices of observation points
% x0_pts = values of height at observation points

[m,n] = size(x0_vals);
I = zeros(m,n);
a = zeros(m,n);
for i=1:max(m,n)
    [a(i), I(i)] = height_val_t(X,X,x0_vals(i));
end
x0_inds = I;
x0_pts = a;

