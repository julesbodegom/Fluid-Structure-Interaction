%% Input data
N   = 20;       % number of cells
Amp = 0;        % amplitude of the initial solution variation
A_m = 0.1;      % amplitude of the mesh motion
u_a = 0;        % advection velocity
tend = 2;       % total simulation time

% Time step refinement parameters
k_vals = 0:6;
dt_array = 0.1 * 2.^(-k_vals);
colors = {'b', 'r', 'g'}; % b=Method 1, r=Method 2, g=Method 3

show_sol = 0; % Keep off for refinement study to speed up execution

% Prepare the log-log plot figure
figure(3);
clf;
hold on;
title('BDF2 Time Step Refinement Study: Error vs \Delta t');
xlabel('\Delta t');
ylabel('Root-Mean-Square Error');
set(gca, 'XScale', 'log', 'YScale', 'log'); % Ensure log-log axes
grid on;

%% Outer loop for the 3 methods
for method = 3
    
    % Arrays to store the computed errors for this method
    err1_array = zeros(1, length(k_vals));
    err2_array = zeros(1, length(k_vals));
    
    %% Inner loop for different time steps
    for k_idx = 1:length(k_vals)
        dt = dt_array(k_idx);
        M = round(tend / dt); % Total number of time steps

        % Initialisation
        dx = zeros(1,N);    % cell volumes
        xi = zeros(1,N+1);  % face centers
        x  = zeros(1,N);    % cell centers

        dx_tn = zeros(1,N);    % cell volumes at tn
        xi_tn = zeros(1,N+1);  % face centers at tn

        dx_tnm1 = zeros(1,N);    % cell volumes at tn-1
        xi_tnm1 = zeros(1,N+1);  % face centers at tn-1

        xi  = (0:N)/N;              % face centers
        xi0 = (0:N)/N;              % face centers at t=0
        x   = (1:N)/N - 0.5/N;      % cell centers
        dx  = xi(2:N+1) - xi(1:N);  % cell volumes

        % determine the initial condition at the cell center locations x
        rho0     = 1 + Amp*sin(x*2*pi); % initial solution
        rho      = rho0;                % set solution to initial solution
        rho_tn   = rho0;                % set solution at tn
        rho_tnm1 = rho0;                % set solution at tn-1

        L = zeros(N);   % Discretization matrix
        
        alpha = [1 -1 0]; % Time discretization coefficients (first step Backward Euler)
        
        sum_sq_err = 0; % Variable to accumulate error for epsilon_2

        %% Simulation loop
        for t=dt:dt:tend
            dx_tn  = dx;  % store cell volume at tn
            xi_tn  = xi;  % store face center at tn
            rho_tn = rho; % store solution at tn
            
            xi     = xi0 + A_m * sin(2*pi*t) * sin(2*pi*xi0); % face centers at tn+1
            dx     = xi(2:N+1)-xi(1:N);                       % cell volumes at tn+1
            x      = xi(1:N) + dx/2;                          % cell centers at tn+1
            
            dxidt_exnp1   = 2*pi*A_m*cos(2*pi*(t))*sin(2*pi*xi0);         % exact face velocity at tn+1
            dxidt_exnp1_2 = 2*pi*A_m*cos(2*pi*(t-dt/2))*sin(2*pi*xi0);    % exact face velocity at tn+1/2
            
            %% IMPLEMENTATION OF DGCL
            if (t == dt)
                dxidt_dgcl = (xi - xi_tn) / dt;  % BE for first step
            else
                dxidt_dgcl = (1.5*xi - 2*xi_tn + 0.5*xi_tnm1) / dt;                                         
            end
            %%

            if (method == 1)
                dxidt = dxidt_exnp1;
            elseif (method ==2)
                dxidt = dxidt_exnp1_2;
            else
                dxidt = dxidt_dgcl;
            end

            % define the relative face velocity
            dxidt_r = u_a - dxidt;

            % Internal cells
            for i=2:N-1
                L(i,i-1:i+1) = [0 alpha(1)*dx(i) 0] + dt * 0.5*[-dxidt_r(i) dxidt_r(i+1)-dxidt_r(i) dxidt_r(i+1)];
            end
            % Boundary cells
            L(1,1:2)   = [alpha(1)*dx(1) 0] + dt * 0.5*[dxidt_r(2)-dxidt_r(1) dxidt_r(2)];
            L(1,N)     =                    - dt * 0.5*dxidt_r(1);
            L(N,N-1:N) = [0 alpha(1)*dx(N)] + dt * 0.5*[-dxidt_r(N) dxidt_r(N+1)-dxidt_r(N)];
            L(N,1)     =                    + dt * 0.5*dxidt_r(1);

            % solve system
            rho = (L\(-alpha(2)*(dx_tn.*rho_tn)'-alpha(3)*(dx_tnm1.*rho_tnm1)'))';

            % update old variables
            dx_tnm1  = dx_tn;
            rho_tnm1 = rho_tn;
            xi_tnm1  = xi_tn;
            
            % After first time step we can use BDF2
            alpha=[3/2 -2 1/2];
            
            % Accumulate error for epsilon_2 (Exact solution is rho0 because Amp=0)
            sum_sq_err = sum_sq_err + sum((rho - rho0).^2);
            
            % Show intermediate solution (disabled for speed)
            if (show_sol)
                figure(1);
                hold off;
                plot(x,rho,'-x');
                axis([0 1 0.8-Amp 1.2+Amp]);
            end
        end

        % Compute final errors for this dt
        err1_array(k_idx) = sqrt( sum((rho - rho0).^2) / N );
        err2_array(k_idx) = sqrt( sum_sq_err / (N * M) );
    end

    %% Plot the errors for this method
    c = colors{method};
    
    % Epsilon 1: Solid line with circle markers
    plot(dt_array, err1_array, [c '-o'], 'DisplayName', ['BDF2 scheme' ' (\epsilon_1)'], 'LineWidth', 1.5);
    
    % Epsilon 2: Dotted line with square markers
    plot(dt_array, err2_array, [c ':s'], 'DisplayName', ['BDF2 scheme' ' (\epsilon_2)'], 'LineWidth', 1.5);
    
end

% Finalize plot
legend('show', 'Location', 'eastoutside');
hold off;