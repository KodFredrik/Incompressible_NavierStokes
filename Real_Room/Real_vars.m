function [t_real, U_real, V_real, W_real, T_real, Xi_real] = Real_vars(t,U,V,W,T,Xi, V0,L,Tmin,Tmax, Ximax)
t_real = t*L/V0;
U_real = U*V0;
V_real = V*V0;
W_real = W*V0;
T_real = T*(Tmax-Tmin) + Tmin;
Xi_real = Xi*Ximax;
end

