

clear all

 addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SWE Bath Data Assimilation 1D'))

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D25/DAB1D25_case_%d.mat',case_num);
load(str)



marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'b--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(case_num)
semilogy([0.01:0.01:0.3],err(end,:),marker1{1},'linewidth',1.5); grid on; hold on;
semilogy([0.01:0.01:0.3],eta_error(end,:),marker2{1},'linewidth',1.5);
set(gca,'FontSize',16);
legend({'$L^2$ Error $\beta^{(b)}$','$L^2$ Error  $\eta(\beta^{(b)})$' },'Interpreter','Latex', 'location', 'southwest');
hold off;
xlabel('Bathymetry Amplitude $\hat{\beta}$','interpreter','latex')
axis square
str2 = sprintf('Bath/Case DAB1D25/Figures/DAB25_case_%d.fig',case_num);
str3 = sprintf('Bath/Case DAB1D25/Figures/DAB25_case_%d.eps',case_num);
savefig(str2);
print(str3, '-depsc2');
end

%%

clear all

 addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SWE Bath Data Assimilation 1D'))

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D25/DAB1D25_case_%d.mat',case_num);
load(str)



marker1 = {'k-+', 'b-o', 'r-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'b--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(1)
semilogy([0.01:0.01:0.3],err(end,:),marker1{case_num},'linewidth',1.5); grid on;hold on;
set(gca,'FontSize',16);

xlabel('Bathymetry Amplitude $\hat{\beta}$','interpreter','latex')
axis square
ylim([1e-3 1e0])
end
hold off;
legend({'Case I','Case II','Case III' },'Interpreter','Latex', 'location', 'northeast');
ylabel('$L^2$ Error $\beta^{(b)}$', 'interpreter','latex')
str2 = 'Bath/Case DAB1D25/Figures/DAB25_bath_err.fig';
str3 = 'Bath/Case DAB1D25/Figures/DAB25_bath_err.eps';
savefig(str2);
print(str3, '-depsc2');

