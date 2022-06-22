% Data Assimilation Scheme for Bathymetry as Control Variable

clear all

for case_num = 1:3

% Fixed parameters for the data assimilation calculations
cfl         = 1/3;
n_obs       = [5 10 20 45];% Number of observation points
rand_x0     = false;        % Use random observation points
ntrial      = 1;            % Number of trials
iter_max    = 3  ;           % Max number of iterations
line_min    = true;        % Find optimal step
nl          = true;         % nonlinear or linear equations
N           = 256 ;          % Number of grid points
xmax        = 3;            % domain size: xmin = -xmax
tmax        = 2*xmax;       % control time
smooth_grad = true;        % Smooth L2 gradient to Sobolev gradient
x0_min      =  0.1*xmax;     % Location of first observation point
%case_num    = 1;
% dmu         = 0.05*xmax;     % Spacing between observation points
% dmu         = abs((xmax-1e-4) - x0_min)/n_obs;
 
%%%%%%%%% Parameters that we are varying %%%%%%
 
amp_bathy = 0.1;            % amplitude of bathymetry (0.3 of full depth is max for periodic ic)

if case_num == 3
    ic_per = true;          % is initial condition periodic (or gaussian)?
    bathy_per = false;           % is bathymetry periodic (or gaussian)?
    tau_n  = -2000;
    filt   = 0.05;          % filtering parameter for Sobolev gradient
    bathy_np_type = 1;          % Gaussian = 1, Sandbar profile =2
elseif case_num == 2
    ic_per = false;          % is initial condition periodic (or gaussian)?
    bathy_per = false;           % is bathymetry periodic (or gaussian)?
    tau_n  = -25000;          % Fixed step
    filt   = 0.05;          % filtering parameter for Sobolev gradient
    bathy_np_type = 2;          % Gaussian = 1, Sandbar profile =2
elseif case_num == 1
    ic_per = false;          % is initial condition periodic (or gaussian)?
    bathy_per = false;           % is bathymetry periodic (or gaussian)?
    tau_n  = -25000;          % Fixed step
    filt   = 0.05;          % filtering parameter for Sobolev gradient  
    bathy_np_type = 1;          % Gaussian = 1, Sandbar profile =2
end
 
amp_ic    = 0.01*amp_bathy; % amplitude of initial condition
k_ic      = 4;              % number of wavelenths of ic in domain
 
k_bathy   = 2;              % number of wavelengths of bathymetry in domain
phi_bathy = pi/2;         % phase shift of bathymetry with respect to initial condition
 
a         = 1;              % range of step size parameter, where the search for next step size is in [tau_n-tau_n/a,tau_n+tau_n/c]
c         =  1/10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conj_grad_type = 0;
iter_chunk     = 10 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisy_obs  = false;  
print_iter = false;        % Print each iteration values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
j_nl=1;
disp('Nonlinear SWE');
    nl=1;
    for k = 1: length(n_obs)
        dmu       = abs((xmax-1e-4) - x0_min)/n_obs(k);
         [eta_exct0(:,k), beta_optimum(:,k), beta_exct0(:,k), X, err(:,k), eta_error(:,k),err_end(k),eta_err_end(k), err_std(:,k), grad(:,k), grad_std(:,k), cost(:,k), cost_std, Y_ex(:,:,k), Y_opt(:,:,k), Y_adj_opt(:,:,k), x0_inds] = ...
           data_assimil_bath_2020_flip(N,n_obs(k),rand_x0,x0_min,dmu, ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, tau_n, amp_bathy, ic_per, bathy_per, ...
    bathy_np_type, amp_ic, k_ic, k_bathy, phi_bathy, conj_grad_type, iter_chunk, a, c, cfl, noisy_obs, print_iter);
    end
    disp(' ');
    
%%
str = sprintf('Bath/Case DAB1D24/DAB1D24_case_%d',case_num);
% Saving Data
save(str, 'ntrial', 'n_obs', 'rand_x0','x0_min', 'dmu', 'line_min','smooth_grad', 'iter_max', ...
    'eta_exct0', 'beta_optimum', 'beta_exct0', 'X', 'err', 'err_std', 'grad', 'eta_error', 'eta_err_end', 'err_end','grad_std', 'cost', 'cost_std', 'Y_ex', 'Y_opt', 'Y_adj_opt', 'x0_inds');

end