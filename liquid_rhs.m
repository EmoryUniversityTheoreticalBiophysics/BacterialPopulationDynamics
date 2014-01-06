function dxdt = liquid_rhs(t,x,K,a,g)

S = x(1);
N = x(2);
dNdt = g*S/(K+S)*N;
dSdt = -dNdt/a;
dxdt = [dSdt; dNdt];
