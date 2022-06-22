clear all
% This script runs the forward SWE equations with:
%  (i) flat bottom
% (ii) batheymetry
% It gives plots at different times of the propagating surface wave

% Fixed parameters for the data assimilation calculations

rand_x0     = false;        % Use random observation points
ntrial      = 1;            % Number of trials
iter_max    = 500;           % Max number of iterations
line_min    = true;        % Find optimal step
tau_n       = 2e3;          % Fixed step
nl          = true;         % nonlinear or linear equations
N           = 1024;          % Number of grid points
xmax        = 3;            % domain size: xmin = -xmax
tmax        = 2*xmax;       % control time
smooth_grad = true;        % Smooth L2 gradient to Sobolev gradient
filt        = 0.05;          % filtering parameter for Sobolev gradient
x0_min      = 0.1*xmax;     % Location of first observation point
% dmu         = 0.05*xmax;     % Spacing between observation points

%%%%%%%%% Parameters that we are varying %%%%%%

amp_bathy = 0.2;            % amplitude of bathymetry (0.3 of full depth is max for periodic ic)
ic_per    = false;          % is initial condition periodic (or gaussian)?
bathy_per = false;           % is bathymetry periodic (or gaussian)?
bathy_np_type = 1;          % Gaussian = 1, Sandbar profile =2
n_obs     = 45;              % Number of observation points
% dmu       = 0.2;
dmu       = abs((xmax-1e-4) - x0_min)/n_obs;
amp_ic    = 0.01*amp_bathy; % amplitude of initial condition
k_ic      = 4;              % number of wavelengths of ic in domain

k_bathy   = 3;              % number of wavelengths of bathymetry in domain
phi_bathy = 0.5*pi;         % phase shift of bathymetry with respect to initial condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Observation points
x0_max = x0_min+(n_obs-1)*dmu;
x0_vals = x0_min:dmu:x0_max;

xmin = -xmax;
h = (xmax - xmin)/N;
sz = h/3;
tmin = 0;

numSteps = round(abs((tmax-tmin)/sz));
tmax = numSteps*sz;
X = xmin:h:xmax-h;
X = reshape(X, [N,1]);

if nl 
disp('Nonlinear SWE')
fw = @fw_SWE_periodic_bath;
bw = @bw_SWE_periodic_mult_bath;
else
disp('Linear SWE')
fw = @fw_SWE_lin_periodic_bath;
bw = @bw_SWE_lin_periodic_mult_bath;
end
    
    
% Exact initial conditions
if (ic_per) 
    disp('Periodic IC')
    eta_exct0 = amp_ic.*cos(k_ic*2*pi/(2*xmax)*X);  
else
    disp('Gaussian IC')
    eta_exct0 = amp_ic.*exp(-X.^2/(0.1*xmax)^2);
end

% Exact Bathymetry
if (bathy_per) 
    disp('Periodic Bathymetry')
    beta_exct0 = amp_bathy*cos(k_bathy*2*pi/(2*xmax)*X+phi_bathy*ones(N,1));
else
     if bathy_np_type == 1
        beta_exct0 = amp_bathy*exp(-(X-1.5).^2/(0.1*xmax)^2);
        disp('Gaussian bathymetry')
     elseif bathy_np_type == 2
        beta_exct0 = (amp_bathy/2)*( tanh(2*(X+2*ones(N,1))) - tanh(2*(X-2*ones(N,1))) );
        disp('Sandbar profile bathymetry')
     else
        disp('Incorrect type entered');
     end
end
disp(['IC amplitude = ', sprintf('%0.7f', amp_ic),', bath amplitude = ', sprintf('%0.3f', amp_bathy)]);


[x0_inds, x0_pts] = mult_x0(X,x0_vals);
% figure(5); plot(X, eta_exct0, X, beta_exct0-1)

% STEP I: Run forward solver using exact bathymetry to get observations
u0 = zeros(N,1);
H = [eta_exct0 ; u0];
[~, ~, obs_x0_mult, ~, T1, Y_ex] = FW_solve_bath(h, sz, H, beta_exct0, 0, tmax, x0_inds, fw);

% STEP iI: Run forward solver using flat bottom case
u0 = zeros(N,1);
H2 = [eta_exct0 ; u0];
beta_exct_zero = zeros(N,1);
[~, ~, obs_x0_mult_zero, ~, T1, Y_ex_zero] = FW_solve_bath(h, sz, H2, beta_exct_zero, 0, tmax, x0_inds, fw);



[Nx,Nt] = size(Y_ex);
time = 1000;
Z = Y_ex(1:Nx/2,time);
Z2 = Y_ex_zero(1:Nx/2,time);

%%

figure (1);
plot(X,Z, 'r', 'linewidth', 1); hold on
plot(X,Z2, 'b','linewidth', 1); hold on
plot(X, 0.01*beta_exct0-0.01, 'k', 'linewidth', 1); 
xlabel('x', 'interpreter', 'latex'); 
set(gca,'FontSize',18);
xlim([min(X) max(X)])
% % ylim([-2e-4, 11e-4])
set(gca,'YTick', [])
hold off;
axis square
print -depsc2 bath/Update_Plots_Feb_2021/eta_bath_no_bath_case_1.eps
savefig('bath/Update_Plots_Feb_2021/eta_bath_no_bath_case_1.fig')


% figure (2);
% % plot(X,Z, 'r', 'linewidth', 1); hold on
% plot(X,abs(Z-Z2), 'k--','linewidth', 1);
% % plot(X, 0.01*beta_exct0-0.01, 'k', 'linewidth', 1); 
% xlabel('x', 'interpreter', 'latex'); 
% set(gca,'FontSize',18);
% xlim([min(X) max(X)])
% % % ylim([-2e-4, 11e-4])
% % set(gca,'YTick', [])
% hold off;
% axis square
% % print -depsc2 bath/err_pert_from_zb_c1.eps


L = length(X);
Y2_a = fft(abs(Z));
P2a = abs(Y2_a/L);
P1a = P2a(1:L/2+1);
P1a(2:end-1) = 2*P1a(2:end-1);
f = (0:(L/2));

Y2_b = fft(abs(Z2));
P2b = abs(Y2_b/L);
P1b = P2b(1:L/2+1);
P1b(2:end-1) = 2*P1b(2:end-1);


figure(2);
plot(f,P1a, 'linewidth', 0.65, 'color', 'b'); hold on; grid on; 
plot(f,P1b, 'linewidth', 0.65, 'color', 'r'); hold on; grid on; 
legend({'Gaussian Bathymetry', 'Flat Bathymetry'},  'interpreter', 'latex','location','northeast');
xlabel('Wavenumber $k$', 'interpreter', 'latex');
ylabel('Energy Spectrum of $\eta(x,t)$', 'interpreter', 'latex');
xlim([0, 30]);
% ylim([0, 1e-4]);
set(gca,'FontSize',18);
hold off
axis square
print -depsc2 bath/Update_Plots_Feb_2021/spec_bath_no_bath_case_1.eps
savefig('bath/Update_Plots_Feb_2021/spec_bath_no_bath_case_1.fig')