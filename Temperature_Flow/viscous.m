function [viscousx, viscousy,viscousz] = viscous(U,V, W, Re,dx,dy,dz)
%VISCOUS Summary of this function goes here
Uxx= diff(U(:,2:end-1, 2:end-1), 2, 1);
Uyy = diff(U(2:end-1,:, 2:end-1), 2, 2);
Uzz = diff(U(2:end-1, 2:end-1, :),2, 3);
viscousx = (1/Re)* (Uxx ./(dx^2) + Uyy./(dy^2) + Uzz./(dz^2));

Vxx= diff(V(:,2:end-1, 2:end-1), 2, 1);
Vyy = diff(V(2:end-1,:, 2:end-1), 2, 2);
Vzz = diff(V(2:end-1, 2:end-1, :),2, 3);
viscousy = (1/Re)* (Vxx./(dx^2) + Vyy./(dy^2) + Vzz./(dz^2));

Wxx= diff(W(:,2:end-1, 2:end-1), 2, 1);
Wyy = diff(W(2:end-1,:, 2:end-1), 2, 2);
Wzz = diff(W(2:end-1, 2:end-1, :),2, 3);
viscousz = (1/Re)* (Wxx./(dx^2) + Wyy./(dy^2) + Wzz./(dz^2));

end

