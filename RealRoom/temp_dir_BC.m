function T = temp_dir_BC(T,X,Y,Z,Size, Tmag)
%This function sets a constant temperature for a square region
sizeT = size(T);

if X == sizeT(1)-1
    horizontal_bounds = Y-floor(Size/2):Y+floor(Size/2);
    vertical_bounds = Z-floor(Size/2):Z+floor(Size/2);
    T(end, horizontal_bounds, vertical_bounds) = Tmag - T(end-1, horizontal_bounds, vertical_bounds);
elseif X == 0
    horizontal_bounds = Y-floor(Size/2):Y+floor(Size/2);
    vertical_bounds = Z-floor(Size/2):Z+floor(Size/2);
    T(1, horizontal_bounds, vertical_bounds) = Tmag - T(2, horizontal_bounds, vertical_bounds);
elseif Y == sizeT(2)-1
    horizontal_bounds = X-floor(Size/2):X+floor(Size/2);
    vertical_bounds = Z-floor(Size/2):Z+floor(Size/2);
    T( horizontal_bounds, end, vertical_bounds) = Tmag - T(horizontal_bounds, end-1, vertical_bounds);
elseif Y == 0
    horizontal_bounds = X-floor(Size/2):X+floor(Size/2);
    vertical_bounds = Z-floor(Size/2):Z+floor(Size/2);
    T( horizontal_bounds, 1, vertical_bounds) = Tmag - T(horizontal_bounds, 2 ,vertical_bounds);
elseif Z == sizeT(3)-1
    horizontal_bounds = X-floor(Size/2):X+floor(Size/2);
    vertical_bounds = Y-floor(Size/2):Y+floor(Size/2);
    T( horizontal_bounds, vertical_bounds, end) = Tmag - T(horizontal_bounds ,vertical_bounds, end-1);
elseif Z == 0
    horizontal_bounds = X-floor(Size/2):X+floor(Size/2);
    vertical_bounds = Y-floor(Size/2):Y+floor(Size/2);
    T( horizontal_bounds, vertical_bounds, 1) = Tmag - T(horizontal_bounds ,vertical_bounds, 2);
end

end

