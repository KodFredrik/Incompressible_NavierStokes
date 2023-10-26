function T = Set_Temperature_BC(T, RoofTemp, FloorTemp)
T = set_scalar_BC(T, 'NNN');
T(:,:,end) = RoofTemp - T(:, end-1 ,:);
T(:,:,1) = FloorTemp - T(:, 2 ,:);
end

