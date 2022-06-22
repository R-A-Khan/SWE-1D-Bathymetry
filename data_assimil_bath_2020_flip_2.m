
function [eta_exct0, beta_optimum, beta_exct0, X, err, eta_error, err_end, eta_err_end, err_std, grad, grad_std, cost, cost_std,  Y_ex, Y_opt, Y_adj_opt, x0_inds] = ...
    data_assimil_bath_2020_flip(N,n_obs,rand_x0,x0_min,dmu,ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, tau_n, amp_bathy, ic_per, bathy_per, ...
    bathy_np_type, amp_ic, k_ic, k_bathy, phi_bathy, conj_grad_type, iter_chunk, a, c, cfl, noisy_obs,  print_iter)

% Input
% % Fixed parameters for the data assimilation calculations
% 
% rand_x0     = false;        % Use random observation points
% ntrial      = 1;            % Number of trials
% iter_max    = 500;           % Max number of iterations
% line_min    = true;        % Find optimal step
% tau_n       = 2e3;          % Fixed step
% nl          = true;         % nonlinear or linear equations
% N           = 128;          % Number of grid points
% xmax        = 3;            % domain size: xmin = -xmax
% tmax        = 2*xmax;       % control time
% smooth_grad = false;        % Smooth L2 gradient to Sobolev gradient
% filt        = 0.1;          % filtering parameter for Sobolev gradient
% x0_min      = 0.1*xmax;     % Location of first observation point
% % dmu         = 0.05*xmax;     % Spacing between observation points
% 
% %%%%%%%%% Parameters that we are varying %%%%%%
% 
% amp_bathy = 0.3;            % amplitude of bathymetry (0.3 of full depth is max for periodic ic)
% ic_per    = false;          % is initial condition periodic (or gaussian)?
% bathy_per = false;           % is bathymetry periodic (or gaussian)?
%                             %(less than or equal 0.3 of depth for periodic initial condition)
% n_obs     = 15;              % Number of observation points
% dmu       = abs((xmax-1e-4) - x0_min)/n_obs
% amp_ic    = 0.01*amp_bathy; % amplitude of initial condition
% k_ic      = 4;              % number of wavelenths of ic in domain
% 
% k_bathy   = 0.5;              % number of wavelengths of bathymetry in domain
% phi_bathy = 0.5*pi;         % phase shift of bathymetry with respect to initial condition
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a         = 2;              % range of step size parameter, where the search for next step size is in [tau_n-tau_n/a,tau_n+tau_n/c]
% c         = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_line_min = false;

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

  for jtrial=1:ntrial
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
    [eta_exact, ~, obs_x0_mult, ~, T1, Y_ex] = FW_solve_bath(h, dt, H, beta_exct0, 0, tmax, x0_inds, fw);
    
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


        % BEGIN LOOP
        iter=0;
        if print_iter
         disp(['Iter','   ','Tau', '   ', '    Cost fminbnd','   ', 'Exit flag','  ', 'cost', '               ',' Grad J','             ','Error', '             ','Eta Error']);
        end
        while norm(grad_J) >= 1e-12 && iter<iter_max
            iter = iter+1;
            % STEP 3
            % Run forward solver to get height at x0 for all t given eta0_1
            H2 = [eta_exct0 ; u0];
            [eta_all, u_all, eta_x0_t_mult, ~, T2, Y1] = FW_solve_bath(h, dt, H2, beta0(:,iter), 0, tmax, x0_inds, fw);

            % STEP 4
            % Run backward solver to get adjoint height at t0 for all x
            H3 = zeros(2*N,1);
            [u_a_t0_x, eta_a_t0_x, ~, Y2] = BW_solve_bath(h, dt, H3, 0, tmax, x0_inds, bw, obs_x0_mult, eta_x0_t_mult, beta0(:,iter), u_all, eta_all, T1, T2);
            % STEP 5: Define gradient of cost function
            u_x_t = Y1(N+1:2*N,:);
            Y2 = fliplr(Y2);
            d_eta_adj_x_t= zeros(N,length(T2));
            for j = 1:length(T2)
                d_eta_adj_x_t(:,j) = cent_diff_u(h,Y2(1:N,j));
            end
 

            prod = u_x_t.*(d_eta_adj_x_t);
            grad_J = trapz(T2,prod,2);
          
            
%            flipY = fliplr(Y2);
%            for j = 1:length(T2)
%              testa(:,j) = cent_diff_u(h,flipY(1:N,j));
%             end
%             grad_J_test = trapz(T2,u_x_t.*testa,2);
%          
            
  
                      


            if smooth_grad
                grad_J = grad_smooth_H2(grad_J,filt);
%                 grad_J_test = grad_smooth_H2(grad_J_test,filt);
            end
            
%           without_flip = norm(grad_J);
%           with_flip = norm(grad_J_test);



            steepest_direction(:,iter) = grad_J;
            % STEP 5.5: Smooth gradient


            conj_direction(:,iter) =  steepest_direction(:,iter); 

            % Plotting J( beta0 - tau.*grad_J(beta0)) as a function of stepsize tau
            tau = -1e5:1000:1e5;
            if plot_line_min && iter == 5
            for i = 1:length(tau)
                J(i) = line_min_mult_bath_CGM(obs_x0_mult, beta0(:,iter), conj_direction(:,iter), N, h, x0_inds, tau(i), eta_exct0, dt, tmin, tmax, fw);
            end
            figure(10);plot(tau, J)
            keyboard
            end
                    
            % STEP 6: Line Minimisation for optimal step size tau_n
            if line_min
            warning('off','optim:fminunc:SwitchingMethod')       
            f = @(tau)line_min_mult_bath_CGM(obs_x0_mult, beta0(:,iter), conj_direction(:,iter), N, h, x0_inds, tau, eta_exct0, dt, tmin, tmax, fw);
%             f = @(tau)line_min_mult_bath_CGM(obs_x0_mult, beta0(:,iter), grad_J_test, N, h, x0_inds, tau, eta_exct0, dt, tmin, tmax, fw);

            % tau_n = fminunc(f,tau0,options);
%             [tau_n, fval, exitflag] = fminbnd(f,tau_n-tau_n/a,tau_n+tau_n/c);
            [tau_n, fval, exitflag] = fminunc(f,tau_n, optimoptions('fminunc','Display','none'));
            %tau_n = min(tau_n, 3/n_obs);
            end


            % STEP 7: Steepest Descent Algorithm
            beta0(:,iter+1)= beta0(:,iter) - tau_n*conj_direction(:,iter);
           
%             beta0(:,iter+1)= beta0(:,iter) + tau_n*(- grad_J_test);
            % beta_optimum = beta0(:,iter+1);
            beta_optimum = beta0(:,iter+1);
            beta_opt_norm = norm(  beta_optimum);
            error = norm(beta_exct0-beta_optimum)/norm(beta_exct0);


            % STEP 8: Run forward solver WITH NEW IC to get height at x0 for all t given
            % (i) optimised Bath   (ii) Exact Bath
            eta_exct2 = amp_ic.*exp(-X.^2/(0.06*xmax)^2);
            H4 = [eta_exct2;u0];
            [eta_opt2, u_opt, eta_x0_t_mult2, ~, T2, Y_opt] = FW_solve_bath(h, dt, H, beta0(:,iter+1), 0, tmax, x0_inds, fw);
            H5 = zeros(2*N,1);
            [~, ~, ~, Y_adj_opt] = BW_solve_bath(h, dt, H5, 0, tmax, x0_inds, bw, obs_x0_mult, eta_x0_t_mult2, beta0(:,iter+1), u_opt, eta_opt2, T1, T2);

                % STEP I: Run forward solver using exact bathymetry to get observations
            u0 = zeros(N,1);
            H_2 = [eta_exct2 ; u0];
            [eta_exact2, ~, ~, ~, ~, ~] = FW_solve_bath(h, dt, H_2, beta_exct0, 0, tmax, x0_inds, fw);
    

            % STEP 8.5: Compute Cost Function
            cost_x = zeros(1,length(T2));
            for i = 1:length(x0_inds)
            cost_x = cost_x + (obs_x0_mult(i,:) - eta_x0_t_mult2(i,:)).^2 ;
            end



            cost(iter,jtrial) = trapz(T2,0.5*cost_x);
            err(iter,jtrial) = norm(beta_exct0-beta_optimum)/norm(beta_exct0);
            grad(iter,jtrial) = norm(grad_J);
            eta_error(iter, jtrial) = sqrt(trapz(T2, trapz( X, (eta_exact2-eta_opt2).^2 ,  1))./trapz(T2, trapz( X, (eta_exact2).^2 ,  1)) ) ;
             
            if print_iter
                if line_min
            disp([sprintf('%1.f',iter),'   ' sprintf('%0.4e',tau_n), '   ' sprintf('%0.7e',fval),'       ' sprintf('%1.0f',exitflag),'          ' sprintf('%0.9e',norm(cost(iter,jtrial))), '          ' sprintf('%0.9e',grad(iter,jtrial)),'          ',sprintf('%0.9e',norm(err(iter,jtrial))),'          ',sprintf('%0.9e',eta_error(iter,jtrial))]);
                else
            disp([sprintf('%1.f',iter),'   ' sprintf('%0.4e',tau_n), '   ' sprintf('%0.7e',norm(cost(iter,jtrial))), '          ' sprintf('%0.7e',grad(iter,jtrial)),'          ',sprintf('%0.7e',norm(err(iter,jtrial)))]);
                end 
            end



                   
        end
  end  


% figure;  
% plot(X, beta_optimum(:,end), X, beta_exct0, 'LineWidth',1.5)
% title(sprintf('N = %d, tau = %d, iterations = %d', N, tau_n, iter_max), 'fontsize',16, 'fontweight', 'bold' )
% set(gca,'FontSize',22)
  
  
err_std  = std(err(end,:));
err      = mean(err,2);

eta_err_std  = std(eta_error(end,:));
eta_error      = mean(eta_error,2);

grad_std = std(grad(end,:));
grad     = mean(grad,2);

cost_std = std(cost(end,:));
cost     = mean(cost,2);

err_end = err(end);
eta_err_end = eta_error(end);


 disp(['N_obs = ', sprintf('%i', n_obs), '  mean error = ', sprintf('%0.4e',err(end)),' mean eta error = ', sprintf('%0.4e',eta_error(end))]);


% 
%  % 
%  figure(6);
% for i = 1:length(T1)
%     % Example of plot
%     Z = Y_ex(1:N,i);
%     Z2 = Y_opt(1:N,i);
%     Z3 = Y_init(1:N,i);
% %     plot(X,Z, 'r', X, Z2, 'b'); %hold on
%     plot(X,Z, 'r'); hold on
%     plot(X, Z2, 'b'); 
%     scatter(X(x0_inds), Y_ex(x0_inds,i), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5); hold off
%     ylim([-2e-3 4e-3])
%     xlim([xmin xmax])
%     title('     Wave propagation error given exact and optimised bathymetry', 'fontsize',8, 'fontweight', 'bold' )
%     set(gca,'FontSize',12)
% %     axis tight
%     % Store the frame
%     M(i)=getframe(gcf);
%     % leaving gcf out crops the frame in the movie.
% end
% 
% 
% % movie(M)
% % fig = figure;
% % movie(gcf,M,1,9)
% 
%   movie2avi(M,'WaveMovie_icnp_bathnp_error.avi')