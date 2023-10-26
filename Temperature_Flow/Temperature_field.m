function T = Temperature_field(T, Ubc,Vbc,Wbc, Re, Pr, dx,dy,dz,dt)

RoofTemp = 0;
FloorTemp = 1;
T_BC = Set_Temperature_BC(T, RoofTemp, FloorTemp);
%Diffusion constant
constant = 1/(Re*Pr);

[TUx, TVy, TWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,T_BC,dx,dy,dz);
Mixterm = TUx(:,2:end-1,2:end-1) + TVy(2:end-1,:,2:end-1) + TWz(2:end-1,2:end-1,:);
DiffT = Diffusion_term(T_BC, constant,dx,dy,dz);
T = T - Mixterm * dt + DiffT * dt;
end

