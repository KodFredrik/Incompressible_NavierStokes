function [t,x,y,z,U,V,W,T,Xi, Q] = Dimensionless_vars ...
    (t_real,x_real,y_real,z_real,U_real,V_real,W_real,T_real,Xi_real,Q_real,...
    V0,L,Tmin,Tmax, Ximax)
%Transforms real quantities to dimensionless form
t = t_real.*V0/L;
x = x_real / L;
y = y_real / L;
z = z_real / L;
U = U_real./V0;
V = V_real./V0;
W = W_real./V0;
T = (T_real - Tmin) ./ (Tmax-Tmin);

Xi = Xi_real ./ Ximax;
Q = (Q_real .* L) ./ (V0*Ximax);
end

