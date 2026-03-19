%% Input data
N   = 20;       % number of cells
Amp = 0;        % amplitude of the initial solution variation
A_m = 0.1;      % amplitude of the mesh motion
u_a = 0;        % advection velocity

dt   = 0.1;     % time step
tend = 2;       % total simulation time

fmt = 'b-';     % linecolor for plot of final solution

show_sol = 1;   % (=1) show intermediate solutions during simulation

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

%% Initialisation

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

% Plot initial solution
figure(1);
hold off;
axis([0 1 0.8-Amp 1.2+Amp]);
plot(x,rho,'-x');
hold off;

L = zeros(N);   % Discretization matrix

alpha = [1 -1 0]; % Time discretization coefficients (first step Backward Euler)

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
  dxidt_exnp1_3 = 2*pi*A_m*cos(2*pi*(t-2*dt))*sin(2*pi*xi0);    % exact face velocity at tn+1/2
  %% IMPLEMENTATION OF DGCL
  if (t == dt)
    dxidt_dgcl = (xi - xi_tn) / dt;  % BE for first step
  else
    dxidt_dgcl    =  (1.5*xi - 2*xi_tn + 0.5*xi_tnm1) / dt;                                         
  end
  %%

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

  % Setting up system to solve; for each cell:
  %
  %    i=3            (rho*dx)^(n+2-i)      
  %   SUM   alpha(i) ----------------- + [ rho_leftface * (dxidt_r)_leftface * (-1) + rho_rightface * (dxidt_r)_rightface * (+1) ]^(n+1) = 0
  %    i=1                 dt                 
  %
  % which can be written for the total system as:
  %
  % L rho^(n+1) = alpha(0) * (rho*dx)^n + alpha(-1) * (rho*dx)^(n-1)
  %
  % We use a simple avarage to compute the solution at a cell face:
  %   rho_face = 0.5 * (rho_leftcell + rho_rightcell)
  %
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
  
  % After first time step we can use 
  alpha=[3/2 -2 1/2];
  
  % Show intermediate solution
  if (show_sol)
    figure(1);
    hold off;
    plot(x,rho,'-x');
    axis([0 1 0.8-Amp 1.2+Amp]);
%    pause
  end
  
end

%% Plot final solution
figure(2);
hold on;
title(['Solution at t=' num2str(tend)]);
ylabel('Solution');
xlabel('x');
plot(x,rho,fmt);
