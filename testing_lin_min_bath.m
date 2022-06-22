clear all
 
% Fixed parameters for the data assimilation calculations
cfl         = 1/3;
n_obs       = [45];% Number of observation points
rand_x0     = false;        % Use random observation points
ntrial      = 1;            % Number of trials
iter_max    = 50  ;           % Max number of iterations
line_min    = true;        % Find optimal step
tau_n       = 8e3;          % Fixed step
nl          = true;         % nonlinear or linear equations
N           = 256 ;          % Number of grid points
xmax        = 3;            % domain size: xmin = -xmax
tmax        = 2*xmax;       % control time
smooth_grad = true;        % Smooth L2 gradient to Sobolev gradient
filt        = 0.05;          % filtering parameter for Sobolev gradient
x0_min      =  0.1*xmax;     % Location of first observation point
% dmu         = 0.05*xmax;     % Spacing between observation points
 dmu       = abs((xmax-1e-4) - x0_min)/n_obs
 
%%%%%%%%% Parameters that we are varying %%%%%%
 
amp_bathy = 0.1;            % amplitude of bathymetry (0.3 of full depth is max for periodic ic)
ic_per    = true;          % is initial condition periodic (or gaussian)?
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




% Observation points
x0_max = x0_min+(n_obs-1)*dmu;
x0_vals = x0_min:dmu:x0_max;

% Spatial domain
xmin = -xmax;
h = (xmax - xmin)/N;
dt = cfl*h; % Time step (CFL = 1)
tmin = 0;

% Assimilation time and step size
numSteps = round(abs((tmax-tmin)/dt));
tmax = numSteps*dt;
X = (xmin:h:xmax-h)';


fw = @fw_SWE_periodic_bath;
bw = @bw_SWE_periodic_mult_bath;

 jtrial=1
     disp(['Trial ', sprintf('%i', jtrial)]);
        if rand_x0
            dx0_min=0;
            while dx0_min < 2*h
                x0_vals = x0_min + (x0_max-x0_min)*rand(n_obs,1);
                for j=1:n_obs-1
                    dx0(j) = x0_vals(j+1)-x0_vals(j);
                end
                dx0_min=min(dx0);
            end
        end

       % Exact initial conditions
    if (ic_per) 
        eta_exct0 = amp_ic.*cos(k_ic*2*pi/(2*xmax)*X);  
        disp('Periodic initial condition')
    else
        eta_exct0 = amp_ic.*exp(-X.^2/(0.1*xmax)^2);
        disp('Gaussian initial condition')
    end

    % Exact Bathymetry
if (bathy_per) 
    disp('Periodic Bathymetry')
    beta_exct0 = amp_bathy*cos(k_bathy*2*pi/(2*xmax)*X+phi_bathy*ones(N,1));
else
     if bathy_np_type == 1
        beta_exct0 = amp_bathy*exp(-(X-0.5*xmax).^2/(0.1*xmax)^2);
        disp('Gaussian bathymetry')
     elseif bathy_np_type == 2
        beta_exct0 = (amp_bathy/2)*( tanh(2*(X+2*ones(N,1))) - tanh(2*(X-2*ones(N,1))) );
        disp('Sandbar profile bathymetry')
     else
        disp('Incorrect type entered');
     end
end

disp(['IC amplitude = ', sprintf('%0.4f', amp_ic),', bath amplitude = ', sprintf('%0.3f', amp_bathy)]);

    % eta_exct0 = eta_exct0 + beta_exct0;
%     plot(X, eta_exct0, X, (beta_exct0-1))
 

    [x0_inds, x0_pts] = mult_x0(X,x0_vals);
    % figure(5); plot(X, eta_exct0, X, beta_exct0-1)
    

    % STEP I: Run forward solver using exact bathymetry to get observations
    u0 = zeros(N,1);
    H = [eta_exct0 ; u0];
    [~, ~, obs_x0_mult, ~, T1, Y_ex] = FW_solve_bath(h, dt, H, beta_exct0, 0, tmax, x0_inds, fw);
    
    % Add gaussian white noise to the observations 
       
%     figure(6);
%     plot(T1,obs_x0_mult(1,:),'r')
%     hold on;
    if noisy_obs
        disp('Noise added to measurements');
        obs_x0_mult = awgn(obs_x0_mult,20,'measured');
    end



    % STEP 2: Distort bathymetry to get "guess"  beta0(:,1) used in algorithm with analytical grad_J
    beta0(:,1) = zeros(N,1);
%     beta0(:,1) = 0.1*beta_exct0;
    
    

    % STEP 2.2: Run forward and backwards solvers to get eta_a_t0_x
    H2 = [eta_exct0;u0];
    [eta_all, u_all, eta_x0_t_mult, ~, T2, Y_init] = FW_solve_bath(h, dt, H2, beta0(:,1), 0, tmax, x0_inds, fw);
    H3 = zeros(2*N,1);
    [u_a_t0_x, eta_a_t0_x, ~, Y2] = BW_solve_bath(h, dt, H3, 0, tmax, x0_inds, bw, obs_x0_mult, eta_x0_t_mult, beta0(:,1), u_all, eta_all, T1, T2);
    Y2 = fliplr(Y2);
    % Define grad_J
    u_x_t = Y_init(N+1:2*N,:);
    d_eta_adj_x_t= zeros(N,length(T2));
    for j = 1:length(T2)
        d_eta_adj_x_t(:,j) = cent_diff_u(h,Y2(1:N,j));
    end

    prod = u_x_t.*d_eta_adj_x_t;
    grad_J = trapz(T2,prod,2);
    norm(grad_J);
    
    if smooth_grad
        disp('Smoothing L2 gradient to H2')
        disp(['filt = ', sprintf('%0.3f', filt)]);
        grad_J = grad_smooth_H2(grad_J,filt);
    end

   

    % Compute Cost Function
    cost_x = zeros(1,length(T2));
    for i = 1:length(x0_inds)
        cost_x = cost_x + (obs_x0_mult(i,:) - eta_x0_t_mult(i,:)).^2 ;
    end

   
    cost(1,jtrial) = trapz(T2,0.5*cost_x) ;
    grad(1,jtrial) = norm(grad_J);
    err(1,jtrial)  = norm(beta_exct0-beta0(:,1))/norm(beta_exct0);


    beta_optimum(:,1) = beta0(:,1);

% Plotting J( beta0 - tau.*grad_J(beta0)) as a function of stepsize tau
tau = -1e5:1000:1e5;

for i = 1:length(tau)
    J(i) = line_min_mult_bath_CGM(obs_x0_mult, beta0(:,1), grad_J, N, h, x0_inds, tau(i), eta_exct0, dt, tmin, tmax, fw);
end
plot(tau, J)