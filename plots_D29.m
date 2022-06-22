clear all

 addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SWE Bath Data Assimilation 1D'))

for case_num = 1
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D29/DAB1D29_case_%d.mat',case_num);
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
str2 = sprintf('Bath/Case DAB1D29/Figures/DAB29_case_%d.fig',case_num);
str3 = sprintf('Bath/Case DAB1D29/Figures/DAB29_case_%d.eps',case_num);
% savefig(str2);
% print(str3, '-depsc2');
end
%%