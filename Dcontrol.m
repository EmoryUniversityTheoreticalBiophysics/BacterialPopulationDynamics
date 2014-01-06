% Numerical solution for non-motile E. coli colony growth as spherical diffusion PDE.
% developed based on ColonyGrowth3D
% N(t) - number of cells in the colony, initially 1 cell
% P(r,t) - concentration of glucose, initially 200ug/ml=2e-10 ug/ml
% D - 50 um^2/s, diffusion coefficient of glucose in agar
% T - endpoint of growth time duration
% R - radius of this spherical model, 1cm = 1e4 * 1um
% g0 - maximum growth rate of E. coli, 1.1 per hour
% a - growth yield, 0.23, a* f0 = dN
% K - Monod constant, constraint for growth rate, 40ug/ml

% Call data instead of calculating codes every time!!!


function [D, DX1, DX2, DX3, DX4, DX5, DX6, DX7] = Dcontrol
T = 35; % hours
R = 6e3; % um
K = 4;  % ug/ml
D0 = 5000; % 50um^2/s = 1.8e5 um^2/h
g0 = 1.1; % per hour
a = 0.92e5; % (ml^-1)/(ug/ml)
dt = 0.001; % hour
dr = 500; % um
nt = T/dt;
nr = R/dr;
V=4/3*pi*R^3;

test=1:1:20;
ntest=length(test);
D=D0*test;
N_Colony = zeros(ntest,nt+1);
N_Liquid = eye(nt+1,1);
Gap = N_Colony;
Rho_Liquid = zeros(1, nt+1);
Rho_Liquid(1) = 500; % ug/ml

for t = 1:nt 
    dN_Liquid = g0*Rho_Liquid(t)/(K+Rho_Liquid(t))*N_Liquid(t)*dt;
    N_Liquid(t+1) = N_Liquid(t)+ dN_Liquid;
    Rho_Liquid(t+1) = Rho_Liquid(t)-dN_Liquid/a ;        
end
Liquidstr = 'Liquid';  
%N_Liquid=N_Liquid*V;
    
    Rho_Colony = zeros(nr+1, nt+1);
    RhoC_r = zeros(nr+1, 1);
    RhoC_rr = zeros(nr+1,1);
    Rho_Colony(:,1) = 500; % ug/ml    
    N_Colony(:,1) = 1;
for i=1:ntest         
          

    for t = 1:nt 
    
        dN_max = g0*Rho_Colony(1,t)/(K+Rho_Colony(1,t))*N_Colony(i,t)*dt; %optimal growth within time dt under P(1,t)
        consumption = dN_max/a ; % resource consumption if growing optimally
    %     S_1=P(1,t)*4/3*pi*dr^3; %local resource supply 
    %     if (consumption > S_1)
    %         dN = S_1*a;
    %         else 
    %         dN = dN_max;
    %     end 
        N_Colony(i,t+1) = N_Colony(i,t) + dN_max;
    
        %calculate dPr
        RhoC_r(1) = (Rho_Colony(2,t)-Rho_Colony(1,t))/dr;
        RhoC_r(2:end-1) = (Rho_Colony(3:end,t)-Rho_Colony(1:end-2,t))/(2*dr);
        RhoC_r(end) = (Rho_Colony(end,t)-Rho_Colony(end-1,t))/dr;
    
        %calculate ddPr
        RhoC_rr(1) = (RhoC_r(2)-RhoC_r(1))/dr;
        RhoC_rr(2:end-1) = (RhoC_r(3:end)-RhoC_r(1:end-2))/(2*dr);
        RhoC_rr(end) = (RhoC_r(end)-RhoC_r(end-1))/dr;
    
        %calculate P(r,t)    
        nflux_1 =  D(i)*RhoC_r(1); %flux from r=2 to 1
        Rho_Colony(1,t+1)=Rho_Colony(1,t)+ (nflux_1*4*pi*dr^2*dt)/(4/3*pi*dr^3)-consumption;
        for r=2:nr+1
            radius=dr*(nr-1);
            diffusion_r_t =  D(i)*(RhoC_r(r)*2/radius + RhoC_rr(r));
            Rho_Colony(r,t+1)=Rho_Colony(r,t)+dt* diffusion_r_t;
        end 
   
    end 
    Gap(i,:) = N_Liquid'-N_Colony(i,:);
    Colonystr{i} = ['Colony: D= ' num2str(D(i)) ' \mum^2/s'];
    gaplabel{i} = ['D = ' num2str(D(i)) ' \mum^2/s'];
end


[Gmax,t_Gmax] = max(Gap,[],2);
tau_Gmax=dt*(t_Gmax-1);
relativegap=zeros(size(Gap));
tau_2=zeros(ntest,1);


for i=1:ntest
    for t=1:nt
        relativegap(i,t)=Gap(i,t)/N_Colony(i,t);        
    end
    tau_2(i)=dt*(find(relativegap(i,:)*100>5, 1)-1);
end 

DX1=N_Liquid; 
DX2=N_Colony; 
DX3=Gap;
DX4=relativegap;
DX7=Gmax;
DX5=tau_Gmax;
DX6=tau_2;


    
    Y = 0:dt:T;
    figure(1)
    semilogy(Y, DX1,Y, DX2,'--');
    xlabel('time (hrs)','FontSize',20)
    ylabel('number of cells','FontSize',20)
    legend([Liquidstr Colonystr],'Location','southeast')
    figure(2)
    plot(Y, DX3);legend(gaplabel,'Location','NorthWest')    
    xlabel('time (hrs)','FontSize',20)
    ylabel('Gap=N_{Liquid}-N_{Colony}','FontSize',20)
    figure(3)
    plot(D,DX6,'o')
    xlabel('D (\mum^2/s)','FontSize',20)
    ylabel('Deviation timepoint(hrs)','FontSize',20)

    
    
