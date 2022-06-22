
clear all
load('resampling.mat','XA','XB','XC')

[N,M] = size(XA);

amp_beta = XA(:,1);
amp_ic   = XA(:,2);
ph_shift = XA(:,3);

tic

for k = 6
    [Y_eta(k), Y_beta(k)] = DA_bath_sobol([amp_beta(k), amp_ic(k), ph_shift(k)]);
    disp(['k = ', sprintf('%i', k),', Bath amplitude = ', sprintf('%0.3f', amp_beta(k)),', IC amplitude = ', sprintf('%0.4f', amp_ic(k)),', Y_beta = ', sprintf('%0.4e', Y_beta(k)),', Y_eta = ', sprintf('%0.4e', Y_eta(k))])
end

toc
% save('YA_sobol.mat', 'Y_eta', 'Y_beta', 'amp_beta', 'amp_ic', 'ph_shift');