function T = Temperature_field(T, Ubc,Vbc,Wbc, constant,windowTemp, radiatorTemp, window_pos, radiator_pos, size, dx,dy,dz,dt)
%This function solves advection-diffusion equation and returns the
%temperature field

%Sets homogeneous neumann condition (adiabatic) everywhere
T_BC= set_scalar_BC(T, 'NNN');

%Set constant dirichlet condition for radiator and window
T_BC  = add_const_temp_BC(T_BC,windowTemp,window_pos,size,'north');
T_BC  = add_const_temp_BC(T_BC,radiatorTemp,radiator_pos,size,'north');

%Computes mixed terms from advective derivative
[TUx, TVy, TWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,T_BC,dx,dy,dz);
%Add them up and exclude boundary points
Mixterm = TUx(:,2:end-1,2:end-1) + TVy(2:end-1,:,2:end-1) + TWz(2:end-1,2:end-1,:);

%Compute laplacian for the diffusion term
DiffT = Diffusion_term(T_BC, constant,dx,dy,dz);

%Step forward in time with explicit Euler
T = T - Mixterm * dt + DiffT * dt;
end

