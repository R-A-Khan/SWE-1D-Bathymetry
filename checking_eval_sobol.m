
load('resampling_ext3.mat')
load('YC3_sobol.mat')

% Find indices of all undefined values
ind = find(isnan(Y_beta));
for i = 1: length(ind)
disp(['k = ', sprintf('%i', ind(i)),', Bath amplitude = ', sprintf('%0.3f', amp_beta(ind(i))),', IC amplitude = ', sprintf('%0.4f', amp_ic(ind(i))),', ph_shift = ', sprintf('%0.4f', ph_shift(ind(i))),', Y_beta = ', sprintf('%0.4e', Y_beta(ind(i))),', Y_eta = ', sprintf('%0.4e', Y_eta(ind(i)))])
end

% Randomly generate indices to remove from Y
% match up sample size to smallest N base rate such that no undefined vals
% in YA, YB, YC

new_length = length(XC3) - length(ind);
for k = 1:2
    ind2(k) = randi(new_length);
end

% Remove indices for undefined values and randomly generated indices
YC3_beta_new = Y_beta;
YC3_beta_new(ind) = [];
YC3_beta_new(ind2) = [];

YC3_eta_new = Y_eta;
YC3_eta_new(ind) = [];
YC3_eta_new(ind2) = [];

XC3_new = XC3;
XC3_new(ind,:) = [];
XC3_new(ind2,:) = [];


% New values should all contain no undefined values, YA and YB should be 
% N x 1, YC should be M*N X 1, corresponding to parameters in
% XA (N x 1), XB (N x 1), XC(M*N x 1)


save('YC3_sobol_new.mat','YC3_beta_new', 'YC3_eta_new', 'XC3_new');