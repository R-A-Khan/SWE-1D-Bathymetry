function J = line_min_mult_bath(obs_x0_mult, beta0, conj_grad, N, h, x0_inds, tau, eta_exct0, dt, tmin, tmax, funsolve)
% Usage:  = line_min_mult(obs_x0_mult, eta0_adj, N, h, x0_inds, tau, eta0, dt, tmin, tmax, funsolve)
%
% Calculates cost function for line minimization algorithm
%
% Input:
% obs_eta_t = observation vector at x0 for all t
% eta0_adj  = adjoint height at t0 for all x
% grad_J    = value assigned to grad_J(beta0)
% tau = stepsize for line minimisation
% eta_exct0 = initial condition at t0 for all x
% dt = time step size
% [tmin, tmax] = time range
% funsolve = function for FD spatial discretisation of forward solver
%
% Output:
% J = cost function

%tau;
N0 = beta0 - tau.*conj_grad;

u = zeros(N,1);
H = [eta_exct0;u];
  [~, ~, eta_x0_t_mult, ~, T2, ~] = FW_solve_bath(h, dt, H, N0, 0, tmax, x0_inds, funsolve);
  


% Define Integrand
sum = zeros(1,length(T2));
for i = 1:length(x0_inds)
    sum = sum + (obs_x0_mult(i,:) - eta_x0_t_mult(i,:)).^2 ;
end
% Numerical Integration using trapezoid Rule
J = trapz(T2,0.5*sum);

% tau_n is argmin of J over tau
% value of tau  where J is a minimum.
end

