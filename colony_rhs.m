function dxdt = colony_rhs(t,x,kappa,a,g,D,j0,dx)

rho = x(1:end-1);
N = x(end);

dNdt = g*rho(j0)/(kappa+rho(j0))*N;

drhodt = zeros(length(rho),1);
% general
drhodt(2:end-1) = D/dx^2*(rho(3:end)+rho(1:end-2)-2*rho(2:end-1));
% endpoints
drhodt(1) = D/dx^2*(rho(2)-rho(1));
drhodt(end) = D/dx^2*(rho(end-1)-rho(end));
% center
drhodt(j0) = drhodt(j0) - 1/dx/a*dNdt;

dxdt = [drhodt; dNdt];
