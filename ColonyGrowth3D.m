% Numerical solution for non-motile E. coli colony growth as spherical diffusion PDE.
% ColonyGrowth3D
% N(t) - number of cells in the colony, initially 1 cell
% P(r,t) - concentration of glucose, initially 200ug/ml
% D - 50 um^2/s, diffusion coefficient of glucose in agar
% T - endpoint of growth time duration
% R - radius of this spherical model, 1cm = 1e4 * 1um
% g0 - maximum growth rate of E. coli, 1.1 per hour
% a - growth yield, 0.23, a* dP(g/L) = dN(g/L)
% K - Monod constant, constraint for growth rate, 4ug/ml
function [N_Colony,N_Liquid,P] = ColonyGrowth3D(P0)
T = 500; % hr
R = 5e3; % um
V=4/3*pi*(R/1e4)^3; % ml
K = 4; % ug/ml
D = 10000; % 50um^2/s = 1.8e4 um^s/hr
g0 = 1.1; % hr^-1
a = 0.92; % (ml^-1)/(ug/ml)
dt = 0.0001; % hr^-1
dr = 500; % um
nt = T/dt;
nr = R/dr;
P = zeros(nr+1, nt+1);
dPr = zeros(nr+1, 1);
ddPr = zeros(nr+1,1);
N_Colony = zeros(1,nt+1);
N_Liquid = zeros(1,nt+1);
Rho_Liquid = zeros(1, nt+1);
P(:,1) = P0; % P0=200ug/ml
Rho_Liquid(1)=P0*V;
N_Colony(1) = 1;
N_Liquid(1)=1;
% f0(nt) = 1/a * g0*P(0,nt)/(K+P(0,nt))*N(nt);

for t = 1:nt 
    dN_Liquid = g0*Rho_Liquid(t)/(K+Rho_Liquid(t))*N_Liquid(t)*dt;
    N_Liquid(t+1) = N_Liquid(t)+ dN_Liquid;
    Rho_Liquid(t+1) = Rho_Liquid(t)-dN_Liquid/(a*2e6) ;        
end
for t = 1:nt 
  if mod(t,nt/100) == 0
    disp([num2str(t/nt*100) '% done'])
  end

    dN_max = g0*P(1,t)/(K+P(1,t))*N_Colony(t)*dt; %optimal growth within time dt under P(1,t)
    consumption = dN_max/a; % resource consumption if growing optimally

    N_Colony(t+1) = N_Colony(t) + dN_max;
    
    %calculate dPr
    dPr(1) = (P(2,t)-P(1,t))/dr;
    dPr(2:end-1) = (P(3:end,t)-P(1:end-2,t))/(2*dr);
    dPr(end) = (P(end,t)-P(end-1,t))/dr;
    
    %calculate ddPr
    ddPr(1) = (dPr(2)-dPr(1))/dr;
    ddPr(2:end-1) = (dPr(3:end)-dPr(1:end-2))/(2*dr);
    ddPr(end) = (dPr(end)-dPr(end-1))/dr;
    
    %calculate P(r,t)    
    nflux_1 =  D*dPr(1); %flux from r=2 to 1
    P(1,t+1)=P(1,t)+ (nflux_1*4*pi*dr^2*dt - consumption)/(4/3*pi*dr^3);
    for r=2:nr+1
        radius=dr*(nr-1);
        diffusion_r_t = dt*D*(dPr(r)*2/radius + ddPr(r));
        P(r,t+1)=P(r,t)+diffusion_r_t;
    end 
   
end 
    X = 0:dr:R;
    Y = 0:dt:T;
    %surf(Y,X,P);
    %rotate3d on
    figure()
    loglog(Y,N_Colony,'--r',Y,N_Liquid,'-g',Y,N_Liquid(end)/100*Y.^(3/2),'--b');
    xlabel('Time(hours)')
    ylabel('Colony size')
% Warning: Axis limits outside float precision, use
% ZBuffer or Painters instead. Not rendering    
% need constraint (P>0)
   
