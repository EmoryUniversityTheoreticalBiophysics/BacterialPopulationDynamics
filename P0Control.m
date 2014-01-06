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


function [X1, X2 ] = P0control
T = 1500; % hours
R = 1e4; % um
K = 4;  % ug/ml
D = 50000; % 50um^2/s = 1.8e5 um^2/h
g0 = 1.1; % per hour
a = 0.92; % (ml^-1)/(ug/ml)
dt = 0.001; % hour
dr = 500; % um
nt = T/dt;
nr = R/dr;
V=4/3*pi*(R/1e4)^3; % ml
Rho_Colony = zeros(nr+1, nt+1);
Rho_Liquid = zeros(1, nt+1);
RhoC_r = zeros(nr+1, 1);
RhoC_rr = zeros(nr+1,1);


test=[1,3,5];
ntest=length(test);
N_Colony = zeros(ntest,nt+1);
N_Liquid = zeros(ntest,nt+1);

for i=1:ntest
Rho_Colony(:,1) = test(i)*100; % ug/ml
Rho_Liquid(1) = test(i)*100; % ug/ml
N_Colony(i,1) = 1;
N_Liquid(i,1) = 1;
% f0(nt) = 1/a * g0*P(0,nt)/(K+P(0,nt))*N(nt);

    for t = 1:nt 
        dN_Liquid = g0*Rho_Liquid(t)/(K+Rho_Liquid(t))*N_Liquid(i,t)*dt;
        N_Liquid(i,t+1) = N_Liquid(i,t)+ dN_Liquid;
        Rho_Liquid(t+1) = Rho_Liquid(t)-dN_Liquid/(a*1e6) ;        
    end
        
    for t = 1:nt 
          if mod(t,nt*ntest/100) == 0
            disp([num2str((t+nt*(i-1))/(nt*ntest)*100) '% done'])
          end
    
        dN_max = g0*Rho_Colony(1,t)/(K+Rho_Colony(1,t))*N_Colony(i,t)*dt; %optimal growth within time dt under P(1,t)
        consumption = dN_max/(a*5e4) ; % resource consumption if growing optimally
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
        nflux_1 =  D*RhoC_r(1); %flux from r=2 to 1
        Rho_Colony(1,t+1)=Rho_Colony(1,t)+ (nflux_1*4*pi*dr^2*dt)/(4/3*pi*dr^3)-consumption;
        for r=2:nr+1
            radius=dr*(nr-1);
            diffusion_r_t =  D*(RhoC_r(r)*2/radius + RhoC_rr(r));
            Rho_Colony(r,t+1)=Rho_Colony(r,t)+dt* diffusion_r_t;
        end 
   
    end 
    Liquidstr{i} = ['Liquid: \rho = ' num2str(Rho_Liquid(1)) ' \mug/ml'];
    Colonystr{i} = ['Colony: \rho = ' num2str(Rho_Colony(1)) ' \mug/ml'];
    gaplabel{i} = ['\rho = ' num2str(Rho_Colony(1)) ' \mug/ml'];
end
%N_Colony = N_Colony/V;
Gap = N_Liquid-N_Colony;
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


X1=N_Liquid; 
X2=N_Colony; 
% X3=Gap;
% X4=relativegap;
% X7=Gmax;
% X5=tau_Gmax;
% X6=tau_2;
% rho=test'*100;


    X = 0:dr:R;
    Y = 0:dt:T;
    figure(1)
	loglog(Y, X1,Y, X2,'--',Y,N_Liquid(end)/100*Y.^(3/2),'--b');
    xlabel('time (hrs)')
    ylabel('number of cells')
    legend([Liquidstr Colonystr],'Location','northeast')
%     figure(2)
%     plot(Y, X3);legend(gaplabel,'Location','NorthWest')
%     figure(3)
%     plot(rho,X6,'o')
%     xlabel('\rho (\mug/ml)')
%     ylabel('Deviation timepoint')
%    
    
