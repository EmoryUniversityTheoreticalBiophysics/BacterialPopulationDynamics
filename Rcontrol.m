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


function [R, RX1, RX2, RX3, RX4, RX5, RX6, RX7] = Rcontrol
T = 35; % hours
%R = 3e4; % um
K = 4;  % ug/ml
D = 50000; % 50um^2/s = 1.8e5 um^2/h
g0 = 1.1; % per hour
a = 0.92; % (ml^-1)/(ug/ml)
dt = 0.001; % hour
dr = 500; % um
nt = T/dt;

test=[1:1:20];
ntest=length(test);
R=3e3*test;
N_Colony = zeros(ntest,nt+1);
N_Liquid = eye(nt+1,1);
Gap = N_Colony;
Rho_Liquid = zeros(1, nt+1);
Rho_Liquid(1) = 500; % ug/ml
    for t = 1:nt 
        dN_Liquid = g0*Rho_Liquid(t)/(K+Rho_Liquid(t))*N_Liquid(t)*dt;
        N_Liquid(t+1) = N_Liquid(t)+ dN_Liquid;
        Rho_Liquid(t+1) = Rho_Liquid(t)-dN_Liquid/(a*1e6) ;        
    end
    Liquidstr = ['Liquid'];    
for i=1:ntest
    
    nr = R(i)/dr;
    Rho_Colony = zeros(nr+1, nt+1);
    RhoC_r = zeros(nr+1, 1);
    RhoC_rr = zeros(nr+1,1);
    Rho_Colony(:,1) = 500; % ug/ml
    
    N_Colony(i,1) = 1;
       

    for t = 1:nt 
    
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
    Gap(i,:) = N_Liquid'-N_Colony(i,:);
    Colonystr{i} = ['Colony: R= ' num2str(R(i)) ' \mum'];
    gaplabel{i} = ['R = ' num2str(R(i)) ' \mum'];
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

RX1=N_Liquid; 
RX2=N_Colony; 
RX3=Gap;
RX4=relativegap;
RX7=Gmax;
RX5=tau_Gmax;
RX6=tau_2;


    
    Y = 0:dt:T;
    figure(1)
	semilogy(Y, RX1,Y, RX2,'--');
    xlabel('time (hrs)','FontSize',20)
    ylabel('number of cells','FontSize',20)
    legend([Liquidstr Colonystr],'Location','southeast')
    figure(2)
    plot(Y, RX3);legend(gaplabel,'Location','NorthWest')    
    xlabel('time (hrs)','FontSize',20)
    ylabel('Gap=N_{Liquid}-N_{Colony}','FontSize',20)
    figure(3)
    plot(R,RX6,'o')
    xlabel('R (\mum)','FontSize',20)
    ylabel('Deviation timepoint','FontSize',20)

    
    
