%% Input data
N   = 20;       % number of cells
Amp = 0;        % amplitude of the initial solution variation
A_m = 0.1;      % amplitude of the mesh motion
u_a = 0;        % advection velocity

tend = 2;       % total simulation time

show_sol = 0;   % (=0) turned off to speed up loop execution

method = 3;     % method = 1 : exact at t(n+1)
                % method = 2 : exact at t(n+1/2)
                % method = 3 : DGCL

if (method == 1)
  disp('Mesh velocities are EXACT at t(n+1)');
elseif (method == 2)
  disp('Mesh velocities are EXACT at t(n+1/2)');
else
  disp('Mesh velocities satisfy DGCL');
end

%% Setup for the 3 requested runs
dt_array  = [0.1, 0.05, 0.025];   % The three time steps
fmt_array = {'b-', 'r-', 'g-'};   % The three plot line styles

% Prepare the final plot figure (Figure 2)
figure(2);
clf; % Clear the figure so you get a fresh plot every time you run the script
hold on;
title(['Solution at t=' num2str(tend) ' for Method ' num2str(method)]);
ylabel('Solution');
xlabel('x');

%% Outer loop for different time steps
for run = 1:length(dt_array)
    dt = dt_array(run);
    fmt = fmt_array{run};
    
    disp(['Running simulation for dt = ', num2str(dt)]);

    %% Initialisation (Reset for each new dt)
    dx   = zeros(1,N);    % cell volumes
    xi   = zeros(1,N+1);  % face centers
    x    = zeros(1,N);    % cell centers

    xi  = (0:N)/N;              % face centers
    xi0 = (0:N)/N;              % face centers at t=0
    x   = (1:N)/N - 0.5/N;      % cell centers
    dx  = xi(2:N+1) - xi(1:N);  % cell volumes

    % determine the initial condition at the cell center locations x
    rho0 = 1 + Amp*sin(x*2*pi); % initial solution
    rho  = rho0;                % set solution to initial solution

    L = zeros(N);   % Discretization matrix

    %% Simulation loop
    for t=dt:dt:tend
      dx_tn  = dx;  % store cell volume at tn
      xi_tn  = xi;  % store face center at tn
      rho_tn = rho; % store solution at tn
      
      xi     = xi0 + A_m * sin(2*pi*t) * sin(2*pi*xi0);   % face centers at tn+1
      dx     = xi(2:N+1)-xi(1:N);                         % cell volumes at tn+1
      x      = xi(1:N) + dx/2;                            % cell centers at tn+1
      
      dxidt_exnp1   = 2*pi*A_m*cos(2*pi*(t))*sin(2*pi*xi0);         % exact face velocity at tn+1
      dxidt_exnp1_2 = 2*pi*A_m*cos(2*pi*(t-dt/2))*sin(2*pi*xi0);    % exact face velocity at tn+1/2
      dxidt_dgcl    = (xi - xi_tn) / dt;                            % face velocity satisfying D-GCL

      if (method == 1)
        dxidt = dxidt_exnp1;
      elseif (method ==2)
        dxidt = dxidt_exnp1_2;
      else
        dxidt = dxidt_dgcl;
      end

      % define the relative face velocity that includes mesh velocity and
      % advection velocity
      dxidt_r = u_a - dxidt;
        
      % Internal cells
      for i=2:N-1
        L(i,i-1:i+1) = [0 dx(i) 0] + dt * 0.5*[-dxidt_r(i) dxidt_r(i+1)-dxidt_r(i) dxidt_r(i+1)];
      end
      % Boundary cells, use periodicity, i.e. index 0 => N, and index N+1 => 1
      L(1,1:2)   = [dx(1) 0] + dt * 0.5*[dxidt_r(2)-dxidt_r(1) dxidt_r(2)];
      L(1,N)     =           - dt * 0.5*dxidt_r(1);
      L(N,N-1:N) = [0 dx(N)] + dt * 0.5*[-dxidt_r(N) dxidt_r(1)-dxidt_r(N)];
      L(N,1)     =           + dt * 0.5*dxidt_r(1);

      % solve system
      rho = (L\(dx_tn.*rho_tn)')';
    end

    %% Plot final solution for this dt
    figure(2);
    plot(x,rho,fmt, 'DisplayName', ['dt = ' num2str(dt)]);
    
end

% Finalize the plot format
figure(2);
legend('show', 'Location', 'best');
hold off;