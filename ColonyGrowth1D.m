% in 1D : Numerical solution for non-motile E. coli colony growth as spherical diffusion PDE.
% ColonyGrowth1D
% N(t) - number of cells in the colony, initially 1 cell
% P(r,t) - concentration of glucose, initially 200ug/ml
% D - 50 um^2/s, diffusion coefficient of glucose in agar
% T - endpoint of growth time duration
% R - radius of this spherical model, 1cm = 1e4 * 1um
% g0 - maximum growth rate of E. coli, 1.1 per hour
% a - growth yield, 0.23, a* dP(g/L) = dN(g/L)
% K - Monod constant, constraint for growth rate, 4ug/ml
function [N,P] = ColonyGrowth1D(P0)
T = 24; % hr
R = 5e3; % um
K = 4; % ug/ml
D = 500; % 50um^2/s = 1.8e4 um^s/hr
g0 = 1.1; % hr^-1
a = 0.92; % (ml^-1)/(ug/ml)
dt = 0.005; % hr^-1
dr = 10; % um
nt = T/dt;
nr = R/dr;
P = zeros(2*nr+1, nt+1);
dPr = zeros(2*nr+1, 1);
ddPr = zeros(2*nr+1,1);
N = zeros(1,nt+1);
P(:,1) = P0; % P0=200ug/ml
N(1) = 1;
% f0(nt) = 1/a * g0*P(0,nt)/(K+P(0,nt))*N(nt);
for t=1:nt 
    dN_max = g0*P(nr+1,t)/(K+P(nr+1,t))*N(t)*dt; %optimal growth within time dt under P(1,t)
    consumption = dN_max/a;
    N(t+1) = N(t) + dN_max;
    
    dPr(1) = (P(2,t)-P(1,t))/dr;
    dPr(2:end-1) = (P(3:end,t)-P(1:end-2,t))/(2*dr);
    dPr(end) = (P(end,t)-P(end-1,t))/dr;
    
    ddPr(1) = (dPr(2)-dPr(1))/dr;
    ddPr(2:end-1) = (dPr(3:end)-dPr(1:end-2))/(2*dr);
    ddPr(end) = (dPr(end)-dPr(end-1))/dr;
    
    P(nr+1,t+1) = P(nr+1,t) + D*ddPr(nr+1)*dt - consumption;
    for r=1:2*nr+1 
        if r==nr+1
            continue
        else P(r,t+1) = P(r,t) + D*ddPr(r)*dt;
        end 
    end 
end 
    X = 0:dr:2*R;
    Y = 0:dt:T;
    %surf(Y,X,P);
    %rotate3d on
    figure(1)
    semilogy(Y,N);
    legend(sprintf('P_0 = %d ug/ml', P0),'Location','NorthWest')
    xlabel('Time(hours)')
    ylabel('Colony size')
    figure(2)
    hold on
    for i=1:nt
        plot(X,P(:,i))
    end 
    hold off
    xlabel('Radius')
    ylabel('Food')

        
