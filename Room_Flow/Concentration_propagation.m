function Xi_new = Concentration_propagation(Xi,constant, Ubc,Vbc,Wbc, Q,dt,dx,dy,dz,pos_out,Size,wall_out)
%CONCENTRATION_PROPAGATION Summary of this function goes here
%   Detailed explanation goes here
XiBC = set_scalar_BC(Xi, 'NNN');
XiBC = Set_Dirichlet_part(XiBC,0,pos_out,Size,wall_out);
%Modify BC for inlet and outlet, let it be zero dirichlet for now.

[XiUx, XiVy, XiWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,XiBC,dx,dy,dz);
Mixterm = XiUx(:,2:end-1,2:end-1) + XiVy(2:end-1,:,2:end-1) + XiWz(2:end-1,2:end-1,:);
DiffXi = Diffusion_term(XiBC, constant,dx,dy,dz);
Xi_new = Xi - Mixterm * dt + DiffXi * dt  + Q * dt;
end

