% Data Assimilation Scheme for Bathymetry as Control Variable
 
clear all
 
% Fixed parameters for the data assimilation calculations
cfl         = 1/3;
n_obs       = [45];% Number of observation points
rand_x0     = false;        % Use random observation points
ntrial      = 1;            % Number of trials
iter_max    = 500  ;           % Max number of iterations
line_min    = true;        % Find optimal step
tau_n       = -25000;          % Fixed step
nl          = true;         % nonlinear or linear equations
N           = 256 ;          % Number of grid points
xmax        = 3;            % domain size: xmin = -xmax
tmax        = 2*xmax;       % control time
smooth_grad = true;        % Smooth L2 gradient to Sobolev gradient
filt        = 0.05;          % filtering parameter for Sobolev gradient
x0_min      =  0.1*xmax;     % Location of first observation point
% dmu         = 0.05*xmax;     % Spacing between observation points
% dmu         = abs((xmax-1e-4) - x0_min)/n_obs;
 
%%%%%%%%% Parameters that we are varying %%%%%%
 
amp_bathy = 0.1;            % amplitude of bathymetry (0.3 of full depth is max for periodic ic)
ic_per    = false;          % is initial condition periodic (or gaussian)?
bathy_per = false;           % is bathymetry periodic (or gaussian)?
                            %(less than or equal 0.3 of depth for periodic initial condition)
bathy_np_type = 1;          % Gaussian = 1, Sandbar profile =2
 
 
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
print_iter = true;        % Print each iteration values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
j_nl=1;
disp('Nonlinear SWE');
    nl=1;
    for k = 1: length(n_obs)
        dmu       = abs((xmax-1e-4) - x0_min)/n_obs(k);
         [eta_exct0(:,k), beta_optimum(:,k), beta_exct0(:,k), X, err(:,k), err_std(:,k), grad(:,k), grad_std(:,k), cost(:,k), cost_std, Y_ex(:,:,k), Y_opt(:,:,k), Y_adj_opt(:,:,k), x0_inds] = ...
           data_assimil_bath_2020_flip(N,n_obs(k),rand_x0,x0_min,dmu, ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, tau_n, amp_bathy, ic_per, bathy_per, ...
    bathy_np_type, amp_ic, k_ic, k_bathy, phi_bathy, conj_grad_type, iter_chunk, a, c, cfl, noisy_obs, print_iter);
    end
    disp(' ');
    
%%

% % Saving Data
% save('Case_DAB1D22iii', 'ntrial', 'n_obs', 'rand_x0','x0_min', 'dmu', 'line_min','smooth_grad', 'iter_max', ...
%     'eta_exct0', 'beta_optimum', 'beta_exct0', 'X', 'err', 'err_std', 'grad', 'grad_std', 'cost', 'cost_std', 'Y_ex', 'Y_opt', 'Y_adj_opt', 'x0_inds');
% 
%% Plotting
% clear all;
% %
% 
% load('Case_DAB1D22iii');

%%
clf;
% iter_max = 80;
iter = 1:iter_max;
marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/20;

set(0,'DefaultFigureWindowStyle','docked')

figure(1);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k)=semilogy(iter(1:skip:end), err(1:skip:end,k), marker1{k},'markersize',10,'linewidth',1.5);hold on;
%          semilogy(iter(1:skip:end), err(1:skip:end,k,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 5', 'N_{obs} = 10', 'N_{obs} = 15', 'N_{obs} = 20', 'N_{obs} = 45'},'location','northeastoutside');
% legend(h, {'N_{obs} = 2','N_{obs} = 3', 'N_{obs} = 4'});
xlabel('Iteration n'); ylabel('$|\beta_{exct}(x) - \beta_{opt}|_2/|\beta_{exct}(x)|_2$','interpreter','latex');xlim([0 iter_max]);
% ylim([1e-2 1e0]);
print -depsc2 err_assimil.eps
hold off;

figure(2);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k)=semilogy(iter(1:skip:end), grad(1:skip:end,k,1)/grad(1,k,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 5', 'N_{obs} = 10', 'N_{obs} = 15', 'N_{obs} = 20', 'N_{obs} = 45'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$|\nabla J^{(n)}(x)|_2/|\nabla J^{(0)}(x)|_2$','interpreter','latex');xlim([0 iter_max]);
print -depsc2 err_grad.eps
hold off;

figure(3);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k) = semilogy(iter(1:skip:end), cost(1:skip:end,k,1)/cost(1,k,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 5', 'N_{obs} = 10', 'N_{obs} = 15', 'N_{obs} = 20', 'N_{obs} = 45'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$J^{(n)}/J^{(0)}$','interpreter','latex');xlim([0 iter_max]);
print -depsc2 err_cost.eps
hold off;
% 
figure(4);
plot(X,beta_optimum(:,1,1),'k-','linewidth',1.5);grid on; hold on;
plot(X,beta_exct0(:,1,1),'k--','linewidth',1.5);grid on; hold on;
set(gca,'FontSize',16);
legend('Optimised','Exact','location','northeastoutside');
xlabel('x'); ylabel('\beta(x)');xlim([min(X) max(X)]);
print -depsc2 optimum_eta.eps
hold off;
% 


figure(5); 
semilogy(X, abs(beta_optimum(:,1,1)-beta_exct0(:,1,1)), 'k-');
set(gca,'FontSize',16); grid on;
xlabel('x'); ylabel('Error'); xlim([min(X) max(X)]);

figure(6);
plot(X,eta_exct0(:,1), 'k-','LineWidth',1.5); hold on
plot(X, 0.1*beta_exct0(:,1,1)-0.1, 'k--','LineWidth',1.5);
scatter(X(x0_inds), Y_ex(x0_inds), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5); hold off
% ylim([-2e-2 4e-3])
xlim([min(X) max(X)])
title('Bathymetry and IC (not to scale)', 'fontsize',8, 'fontweight', 'bold' )
set(gca,'FontSize',12)
    
% [Nx,Nt] = size(Y_opt);
%  figure(6);
% for i = 1:10:Nt
%     % Example of plot
%     Z = Y_ex(1:Nx/2,i);
% %     Z2 = Y_opt(1:Nx/2,i);
% %     plot(X,Z, 'r', X, Z2, 'b'); %hold on
%     plot(X,Z, 'r'); hold on
%     plot(X, 0.01*beta_exct0-0.01, 'b'); 
%     scatter(X(x0_inds), Y_ex(x0_inds,i), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5); hold off
%     ylim([-4e-2 4e-3])
%    xlim([min(X) max(X)])
%     title('     Wave propagation error given exact and optimised bathymetry', 'fontsize',8, 'fontweight', 'bold' )
%     set(gca,'FontSize',12)
% %     axis tight
%     % Store the frame
%     M(i)=getframe(gcf);
%     % leaving gcf out crops the frame in the movie.
% end