function Xi_new = Concentration_propagation(Xi,constant, U,V,W, Q,dt,dx,dy,dz)
%CONCENTRATION_PROPAGATION Summary of this function goes here
%   Detailed explanation goes here
XiBC = set_scalar_BC(Xi, 'NNN');
[U,V,W] = set_Dirichlet_BC(U,V,W);
[XiUx, XiVy, XiWz] = Scalar_mixed_terms(U,V,W,XiBC,dx,dy,dz);
Mixterm = XiUx(:,2:end-1,2:end-1) + XiVy(2:end-1,:,2:end-1) + XiWz(2:end-1,2:end-1,:);
DiffXi = Diffusion_term(XiBC, constant,dx,dy,dz);
Xi_new = Xi - Mixterm * dt + DiffXi * dt  + Q * dt;
end

