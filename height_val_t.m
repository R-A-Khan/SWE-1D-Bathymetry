function [a,I] = height_val_t(f, t_vec, t_val)
% Usage:
% [a,I] = height_val_t(f, t_vec, t_val)
%
% Take any time value as input, and compute the corresponding height at
% that time (or closest value to it in temporal discretisation)
%
% Input:
% f = height vector at all t
% t_val = value of time input
% t_vec = discrete time steps corresponding to f
%
% Output:
% a = height values
% I = positions of observations


% Find index of t_vec value closest to t_val
tmp = abs(t_vec-t_val); 
[~, I] = min(tmp);

a = f(I) ;

