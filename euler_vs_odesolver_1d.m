clear all

% parameters
N0 = 1;
S0 = 100; % ug
a = 1e6; % 1/ug
g = .7; % 1/hr
T1 = 1e3; % hr
T2 = 1e5; % hr
K = S0/10; % ug
D = .1; % um^2/s
X = 50; % um
kappa = K/(2*X); % ug/um



% 1. liquid: Euler method (fixed time step)

% define time step
dt = (1/g)/1e3; % hr
Nt = round(T1/dt);
t1 = 0:dt:T1;

% initialize
N1 = [N0; zeros(Nt,1)];
S1 = [S0; zeros(Nt,1)];

% integrate
rt0 = cputime;
for i = 1:Nt
  N1(i+1) = N1(i)*(1+dt*g*S1(i)/(K+S1(i)));
  S1(i+1) = S1(i) - dt*(1/a)*g*S1(i)/(K+S1(i))*N1(i);
end
rt1 = cputime-rt0;

% for legend
lstr{1} = sprintf('liquid: Euler (runtime = %0.1f s)',rt1);





% 2. liquid: Matlab ode solver (adaptive time step)

% initialize vector x = [S; N]
x0 = [S0; N0];

% integrate
rt0 = cputime;
[t2,x] = ode15s(@(t,x)liquid_rhs(t,x,K,a,g),[0 T2],x0);
rt2 = cputime-rt0;

% extract N(t) vector
N2 = x(:,2);

% for legend
lstr{2} = sprintf('liquid: ode15s (runtime = %0.1f s)',rt2);





% 3. colony: Euler method (fixed time step)

% define discretization in space
dx = X/30; % um
Nx = round(2*X/dx);
j0 = ceil(Nx/2); % index of position of colony

% initialize
t3 = 0:dt:T1;
N3 = [N0; zeros(Nt,1)];
rho0 = S0/(2*X)*ones(1,Nx); % ug/um (uniform)
rho3 = [rho0; zeros(Nt,Nx)];

% integrate
rt0 = cputime;
for i = 1:Nt
  N3(i+1) = N3(i)*(1+dt*g*rho3(i,j0)/(kappa+rho3(i,j0)));
  rho3(i+1,2:end-1) = rho3(i,2:end-1) + D/(dx^2/dt)...
      *(rho3(i,3:end)+rho3(i,1:end-2)-2*rho3(i,2:end-1));
  % endpoints
  rho3(i+1,1) = rho3(i,1) ...
      + D/(dx^2/dt)*(rho3(i,2)-rho3(i,1));
  rho3(i+1,end) = rho3(i,end) ...
      + D/(dx^2/dt)*(rho3(i,end-1)-rho3(i,end));
  % center
  rho3(i+1,j0) = rho3(i+1,j0) ...
      - dt/dx*(1/a)*g*rho3(i,j0)/(kappa+rho3(i,j0))*N3(i);
end
rt3 = cputime-rt0;

% for legend
lstr{3} = sprintf('colony: Euler (runtime = %0.1f s)',rt3);





% 4. colony: Matlab ode solver (adaptive time step)

% initialize vector x = [rho; N]
x0 = [rho0'; N0];

% integrate
rt0 = cputime;
[t4,x] = ode15s(@(t,x)colony_rhs(t,x,kappa,a,g,D,j0,dx),[0 T2],x0);
rt4 = cputime-rt0;

% extract N(t) vector
N4 = x(:,end);

% for legend
lstr{4} = sprintf('colony: ode15s (runtime = %0.1f s)',rt4);





% scaling
t_ = [1e2 1e4];
N_ = 10^5.5*t_.^.5;
lstr{5} = '\propto t^{1/2}';


% plot
figure(1); clf
loglog(t1,N1,'b-',t2,N2,'r--',t3,N3,'g-',t4,N4,'m--',t_,N_,'k-')
xlim([1e-1 T2])
ylim([1e0 1e9])
xlabel('time (hours)')
ylabel('number of cells')
legend(lstr,'location','se')

