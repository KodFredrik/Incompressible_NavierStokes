function [nonlinU, nonlinV, nonlinW] = nonlinear(U,V,W, dx,dy,dz)
%This function computes all the non-linear terms for the velocities
%It is done by averaging the quantities to coincide at the control surface
%where the quantity is being calculated. The quantities are then multiplied
%and differentiated. The function then returns the values for the
%non-linear terms at the interior of the domain.

%Calculate non-linear terms, begin with squared terms
Uavgx = average(U,1);
Usquared = Uavgx.^2;

Vavgy = average(V,2);
Vsquared = Vavgy.^2;

Wavgz = average(W,3);
Wsquared = Wavgz.^2;

%mixed term UV
Uavgy = average(U,2);
Vavgx = average(V,1);
UV = Uavgy.*Vavgx;

%mixed term UW
Uavgz = average(U,3);
Wavgx = average(W,1);
UW = Uavgz.*Wavgx;

%mixed term VW
Vavgz = average(V,3);
Wavgy = average(W,2);
VW = Vavgz .* Wavgy;

%Derivatives for nonlinear terms for U
USQRDx = (1/dx) * diff(Usquared,1,1);
UVy =(1/dy) * diff(UV,1,2);
UWz = (1/dz) * diff(UW,1,3);

%Only include interior points (2:end-1 to exclude outside points in
%directions where derivative has not been taken
nonlinU = USQRDx(:,2:end-1,2:end-1) + UVy(2:end-1,:,2:end-1) + UWz(2:end-1,2:end-1,:);

%Derivatives for nonlinear terms for V
VSQRDy = (1/dy) * diff(Vsquared,1,2);
UVx = (1/dx) * diff(UV,1,1);
VWz = (1/dz) * diff(VW,1,3);

nonlinV = VSQRDy(2:end-1,:,2:end-1) + UVx(:,2:end-1,2:end-1) + VWz(2:end-1,2:end-1,:);

%Derivatives for the nonlinear terms for W
WSQRDw = (1/dz) * diff(Wsquared,1,3);
UWx = (1/dx) * diff(UW,1,1);
VWy = (1/dy) * diff(VW,1,2);

nonlinW = WSQRDw(2:end-1,2:end-1,:) + UWx(:,2:end-1,2:end-1) + VWy(2:end-1, :, 2:end-1);
end

