function Xi_new = Concentration_propagation(Xi,constant, Ubc,Vbc,Wbc, Q,dt,dx,dy,dz,pos_out,Size_out,wall_out, pos_in, Size_in, wall_in)
%CONCENTRATION_PROPAGATION Summary of this function goes here
%   Detailed explanation goes here
XiBC = set_scalar_BC(Xi, 'NNN');
%Set outlet inlet to zero dirichlet to act as sink
XiBC = Set_Dirichlet_part(XiBC,0,pos_out,Size_out,wall_out);
XiBC = Set_Dirichlet_part(XiBC,0,pos_in,Size_in,wall_in);
%Size of neumann condition? V0*size*dt multiply by conc??

%XiBC = set_scalar_neumann_part(XiBC,dx^3,pos_out,Size_out,wall_out,dx,dy,dz);
%XiBC = set_scalar_neumann_part(XiBC,dx^3,pos_in,Size_in,wall_in,dx,dy,dz);


[XiUx, XiVy, XiWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,XiBC,dx,dy,dz);
Mixterm = XiUx(:,2:end-1,2:end-1) + XiVy(2:end-1,:,2:end-1) + XiWz(2:end-1,2:end-1,:);
DiffXi = Diffusion_term(XiBC, constant,dx,dy,dz);
Xi_new = Xi - Mixterm * dt + DiffXi * dt  + Q * dt;
end

