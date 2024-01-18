function Xi_new = Concentration_propagation(Xi,constant, Ubc,Vbc,Wbc, Q,dt,dx,dy,dz,pos_out,Size_out,wall_out, pos_in, Size_in, wall_in)
%This function solves the advection-diffusion equation for the propagation
%of a concentration.

%Sets homogenous neumann conditions (adiabatic) everywhere
XiBC = set_scalar_BC(Xi, 'NNN');

%Set outlet inlet to zero dirichlet to act as sink
XiBC = Set_Dirichlet_part(XiBC,0,pos_out,Size_out,wall_out);
XiBC = Set_Dirichlet_part(XiBC,0,pos_in,Size_in,wall_in);

%Computes mixed terms from advective derivative
[XiUx, XiVy, XiWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,XiBC,dx,dy,dz);
%Add them up and exclude boundary points
Mixterm = XiUx(:,2:end-1,2:end-1) + XiVy(2:end-1,:,2:end-1) + XiWz(2:end-1,2:end-1,:);

%Compute laplacian for the diffusion term
DiffXi = Diffusion_term(XiBC, constant,dx,dy,dz);

%Step forward in time with explicit Euler
Xi_new = Xi - Mixterm * dt + DiffXi * dt  + Q * dt;
end

