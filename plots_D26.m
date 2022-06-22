

clear all

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D26/DAB1D26_case_%d.mat',case_num);
load(str)

if case_num == 1;
err(:,11:14) = [];
eta_error(:,11:14) = [];
end

marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'b--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(case_num)
semilogy([0.0001:0.0001:0.001],err(end,:),marker1{1},'linewidth',1.5); grid on; hold on;
semilogy([0.0001:0.0001:0.001],eta_error(end,:),marker2{1},'linewidth',1.5);
set(gca,'FontSize',16);
legend({'$L^2$ Error $\beta^{(b)}$','$L^2$ Error  $\eta(\beta^{(b)})$' },'Interpreter','Latex', 'location', 'northeast');
hold off;
xlabel('$\hat{\eta}$','interpreter','latex')
axis square
str2 = sprintf('Bath/Case DAB1D26/Figures/DAB26_case_%d.fig',case_num);
str3 = sprintf('Bath/Case DAB1D26/Figures/DAB26_case_%d.eps',case_num);
savefig(str2);
print(str3, '-depsc2');

end
