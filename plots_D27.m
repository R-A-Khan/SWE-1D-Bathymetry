close all 
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results with smoothing for Kappa test

load('bath/Case DAB1D9/Case_DAB1D9i');
kappa_filt(:,1)= kappa;
load('bath/Case DAB1D9/Case_DAB1D9iv');
kappa_filt(:,3)= kappa;
load('bath/Case DAB1D9/Case_DAB1D9v');
kappa_filt(:,2)= kappa;

marker1 = {'k-+', 'b-o', 'r-x', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'b--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip = 5

set(0,'DefaultFigureWindowStyle','docked')
ep = logspace(-8,-1,100);
N = 512;
filt = 0.05;
l = 1/( round(filt*(N/2)));

for l = 1:3
figure(1);
loglog(ep(1:skip:length(kappa_filt(:,l))), abs(abs(kappa_filt(1:skip:end,l))-1),marker1{l}, 'linewidth',1.5); grid on; hold on
set(gca,'FontSize',20)
% ylim([(8e-3), (1.2e-2)])
%title('Kappa test for Numerical  Solution  II')
end
xlabel('$log(\epsilon)$', 'interpreter','latex');
ylabel('$|\kappa(\epsilon)-1|$','interpreter','latex');
legend({'Case I','Case II','Case III' },'Interpreter','Latex', 'location', 'northeast');
hold off;
axis square
str2 = 'Bath/Case DAB1D27/DAB27_kappa.fig';
str3 = 'Bath/Case DAB1D27/DAB27_kappa.eps';
savefig(str2);
print(str3, '-depsc2');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results with smoothing for DA

marker1 = {'k-+', 'b--', 'r-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D27/DAB1D27_case_%d.mat',case_num);
load(str)

iter = 1:iter_max;
skip=iter_max/20;


figure(2);
% semilogy(iter(1:skip:end), err_no_filt_iv(1:skip:end,1), marker1{1},'markersize',10,'linewidth',1.5); grid on; hold on;
semilogy(iter(1:skip:end), err(1:skip:end), marker1{case_num},'markersize',10,'linewidth',1.5);hold on;
% semilogy(iter(1:skip:end), err_filt(1:skip:end,1), marker2{1},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',18); grid on;
xlabel('Iteration n', 'interpreter', 'latex'); ylabel('Relative $L^2$ error for $\beta(x)$','interpreter','latex');
xlim([0 iter_max]);
legend({'Case I','Case II','Case III' },'Interpreter','Latex', 'location', 'northeast');
% ylim([1e-2 1e0]);
hold off;
axis square
str2 = 'Bath/Case DAB1D27/DAB27_err.fig';
str3 = 'Bath/Case DAB1D27/DAB27_err.eps';
savefig(str2);
print(str3, '-depsc2');


%%

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D27/DAB1D27_case_%d.mat',case_num);
load(str)

iter = 1:iter_max;
skip=iter_max/20;
marker1 = {'k-+', 'b--', 'r-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(3);
semilogy(iter(1:skip:end), cost(1:skip:end)/cost(1), marker1{case_num},'markersize',10,'linewidth',1.5);hold on; grid on;
end
legend({'Case I','Case II','Case III' },'Interpreter','Latex', 'location', 'northeast');
xlabel('x', 'interpreter', 'latex');ylabel('$J^{(n)}/J^{(0)}$','interpreter','latex');xlim([0 iter_max]);
set(gca,'FontSize',18);
hold off;
axis square

%%
str2 = 'Bath/Case DAB1D27/DAB27_cost.fig';
str3 = 'Bath/Case DAB1D27/DAB27_cost.eps';
savefig(str2);
print(str3, '-depsc2');

%%

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D27/DAB1D27_case_%d.mat',case_num);
load(str)

figure(case_num);
plot(X,beta_optimum,'k-','linewidth',1.5);grid on; hold on;
plot(X,beta_exct0,'k--','linewidth',1.5, 'color', 'r');grid on;
% plot(X,beta_opt_filt(:,1,1),'k--','linewidth',1.5);grid on; 

set(gca,'FontSize',18);
% legend('Optimised No Smoothing', 'Optimised with Smoothing','location','northeast');
xlabel('x', 'interpreter', 'latex'); ylabel('$\beta(x)$', 'interpreter', 'latex');xlim([min(X) max(X)]);
legend({'$H^2$ Smoothing','Exact'},'interpreter', 'latex','location','northwest');
% print -depsc2 optimum_eta.eps

hold off;
axis square
str2 = sprintf('Bath/Case DAB1D27/DAB27_opt_bath_%d.fig', case_num);
str3 = sprintf('Bath/Case DAB1D27/DAB27_opt_bath_%d.eps', case_num);
savefig(str2);
print(str3, '-depsc2');
end


