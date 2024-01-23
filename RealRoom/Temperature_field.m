function T = Temperature_field(T, Ubc,Vbc,Wbc, constant, dx,dy,dz,dt, temp_dir_array)
%This function solves advection-diffusion equation and returns the
%temperature field

%Sets homogeneous neumann condition (adiabatic) everywhere
T_BC= set_scalar_BC(T, 'NNN');

%Set constant dirichlet condition for radiator and window
tempdir_size = size(temp_dir_array);
for i = 1:tempdir_size(1)
    T_BC = temp_dir_BC(T_BC, temp_dir_array(i,1),temp_dir_array(i,2),temp_dir_array(i,3),temp_dir_array(i,4),temp_dir_array(i,5));
end
%T_BC  = add_const_temp_BC(T_BC,windowTemp,window_pos,size,window_wall);
%T_BC  = add_const_temp_BC(T_BC,radiatorTemp,radiator_pos,size,radiator_wall);

%Computes mixed terms from advective derivative
[TUx, TVy, TWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,T_BC,dx,dy,dz);
%Add them up and exclude boundary points
Mixterm = TUx(:,2:end-1,2:end-1) + TVy(2:end-1,:,2:end-1) + TWz(2:end-1,2:end-1,:);

%Compute laplacian for the diffusion term
DiffT = Diffusion_term(T_BC, constant,dx,dy,dz);

%Step forward in time with explicit Euler
T = T - Mixterm * dt + DiffT * dt;
end

