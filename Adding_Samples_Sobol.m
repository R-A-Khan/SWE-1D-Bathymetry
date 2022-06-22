
load('YA_sobol_new.mat')
load('YB_sobol_new.mat')
load('YA2_sobol_new.mat')
load('YB2_sobol_new.mat')
XA12 = vertcat(XA_new, XA2_new);
XB12 = vertcat(XB_new, XB2_new);

N = length(XA12);
N2 = 3000 ; % increase of base sample size
% (that means: N2*(M+2) new samples that will need to be evaluated)
X(1:N,:) = XA12;
X(N+1:2*N,:) = XB12;
Xext2 = AAT_sampling_extend(X,DistrFun,DistrPar,2*(N+N2)) ; % extended sample 
% (it includes the already evaluated samples 'X' and the new ones
Xnew = Xext2(2*N+1:end,:) ; % extract the new input samples that need to be evaluated
% Resampling strategy:
[ XA3, XB3, XC3 ] = vbsa_resampling(Xnew);
save('resampling_ext3.mat','XA3','XB3','XC3')