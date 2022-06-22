


set(0,'DefaultFigureWindowStyle','docked')
case_num = 3;
str = sprintf('Bath/Case DAB1D24/DAB1D24_case_%d.mat',case_num);

load(str)
for k = 1:4
    
iter = 1:iter_max;
marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/20;

figure(k);
semilogy(iter(1:skip:end),err(1:skip:end,k), marker1{k},'linewidth',1.5);grid on; hold on;
semilogy(iter(1:skip:end),eta_error(1:skip:end,k),marker2{k},'linewidth',1.5);
xlim([0 iter_max]);
set(gca,'FontSize',16);
legend({'$L^2$ Error $\beta^{(b)}$','$L^2$ Error  $\eta(\beta^{(b)})$' },'Interpreter','Latex', 'location', 'northeast');
hold off;
axis square

end

%%
clear all

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D24/DAB1D24_case_%d.mat',case_num);
load(str)


marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'b--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(case_num)
semilogy([5 10 20 45],err(end,:),marker1{1},'linewidth',1.5); grid on; hold on;
semilogy([5 10 20 45],eta_error(end,:),marker2{1},'linewidth',1.5);
set(gca,'FontSize',16);
legend({'$L^2$ Error $\beta^{(b)}$','$L^2$ Error  $\eta(\beta^{(b)})$' },'Interpreter','Latex', 'location', 'northeast');
hold off;
xlabel('$N_{obs}$','interpreter','latex')
axis square
str2 = sprintf('Bath/Case DAB1D24/Figures/DAB24_case_%d.fig',case_num);
str3 = sprintf('Bath/Case DAB1D24/Figures/DAB24_case_%d.eps',case_num);
savefig(str2);
print(str3, '-depsc2');
end


%% Number of observation points

set(0,'DefaultFigureWindowStyle','docked')
for case_num = 1:3;
str = sprintf('Bath/Case DAB1D24/DAB1D24_case_%d.mat',case_num);

load(str)
for k = 1:4
    
iter = 1:iter_max;
marker1 = {'k-+', 'b-o', 'r-*', 'g-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/20;

figure(case_num);
semilogy(iter(1:skip:end),err(1:skip:end,k), marker1{k},'linewidth',1.5);grid on; hold on;
xlim([0 iter_max]);
set(gca,'FontSize',16);
legend({'$N_{obs} = 5$','$N_{obs} = 10$','$N_{obs} = 20$','$N_{obs} = 45$' },'Interpreter','Latex', 'location', 'northeast');
xlabel('Iteration n', 'interpreter', 'latex'); ylabel('Relative $L^2$ error for $\beta(x)$','interpreter','latex');
axis square
end

str3 = sprintf('Bath/Case DAB1D24/Figures/D24_obs_%d_err.eps', case_num);
print(str3, '-depsc2')
str5 = sprintf('Bath/Case DAB1D24/Figures/D24_obs_%d_err.fig', case_num);
savefig(str5)
end


%%

set(0,'DefaultFigureWindowStyle','docked')
for case_num = 1:3;
str = sprintf('Bath/Case DAB1D24/DAB1D24_case_%d.mat',case_num);

load(str)
for k = 1:4
    
iter = 1:iter_max;
marker1 = {'k-+', 'b-o', 'r-*', 'g-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/20;

figure(case_num);
semilogy(iter(1:skip:end),cost(1:skip:end,k)/cost(1,k), marker1{k},'linewidth',1.5);grid on; hold on;
xlim([0 iter_max]);
set(gca,'FontSize',16);
legend({'$N_{obs} = 5$','$N_{obs} = 10$','$N_{obs} = 20$','$N_{obs} = 45$' },'Interpreter','Latex', 'location', 'northeast');
xlabel('Iteration n', 'interpreter', 'latex');
ylabel('$J^{(n)}/ J^{(0)}$','interpreter','latex');
axis square
end

str3 = sprintf('Bath/Case DAB1D24/Figures/D24_obs_%d_cost.eps', case_num);
print(str3, '-depsc2')
str5 = sprintf('Bath/Case DAB1D24/Figures/D24_obs_%d_cost.fig', case_num);
savefig(str5)
end

%%

for case_num = 1:3
set(0,'DefaultFigureWindowStyle','docked')
str = sprintf('Bath/Case DAB1D24/DAB1D24_case_%d.mat',case_num);

load(str)
marker1 = {'k-', 'b-o', 'r-*', 'g-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'r--','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(case_num);
plot(X,beta_optimum(:,1), marker1{1},'linewidth',1.5);grid on; hold on;
plot(X,beta_optimum(:,4), marker2{1},'linewidth',1.5);grid on;
% plot(X,beta_opt_filt(:,1,1),'k--','linewidth',1.5);grid on; 

set(gca,'FontSize',18);
% legend('Optimised No Smoothing', 'Optimised with Smoothing','location','northeast');
xlabel('x', 'interpreter', 'latex'); ylabel('$\beta(x)$', 'interpreter', 'latex');xlim([min(X) max(X)]);
legend({'$N_{obs} = 5$','$N_{obs} = 45$'},'interpreter', 'latex','location','northwest');
% print -depsc2 optimum_eta.eps

hold off;
axis square
str2 = sprintf('Bath/Case DAB1D24/DAB24_opt_bath_N_obs_%d.fig', case_num);
str3 = sprintf('Bath/Case DAB1D24/DAB24_opt_bath_N_obs_%d.eps', case_num);
savefig(str2);
print(str3, '-depsc2');
end



