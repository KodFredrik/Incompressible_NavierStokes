function T = Temperature_field(T, Ubc,Vbc,Wbc, constant,windowTemp, radiatorTemp, window_pos, radiator_pos, size, dx,dy,dz,dt)
T_BC= set_scalar_BC(T, 'NNN');
%T_BC(:,:,end) = RoofTemp - T_BC(:, : ,end-1);
%T_BC(:,:,1) = FloorTemp - T_BC(: ,:, 2);

%T_BC(end,:,:) = RoofTemp - T_BC(end-1, : ,:);
%T_BC(1,:,:) = FloorTemp - T_BC(2 ,:, :);

%Add cold window over hot radiator on northern wall

T_BC  = add_const_temp_BC(T_BC,windowTemp,window_pos,size,'north');
T_BC  = add_const_temp_BC(T_BC,radiatorTemp,radiator_pos,size,'north');

[TUx, TVy, TWz] = Scalar_mixed_terms(Ubc,Vbc,Wbc,T_BC,dx,dy,dz);
Mixterm = TUx(:,2:end-1,2:end-1) + TVy(2:end-1,:,2:end-1) + TWz(2:end-1,2:end-1,:);
DiffT = Diffusion_term(T_BC, constant,dx,dy,dz);
T = T - Mixterm * dt + DiffT * dt;
end

