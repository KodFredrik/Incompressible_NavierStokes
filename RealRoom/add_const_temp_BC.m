function X = add_const_temp_BC(X,Xmag,pos,Size,wall)
%This function sets a constant temperature for a square region

horizontal_bounds = pos(1)-floor(Size/2):pos(1)+floor(Size/2);
vertical_bounds = pos(2)-floor(Size/2):pos(2)+floor(Size/2);

if strcmp(wall,'east')
    X(end, horizontal_bounds, vertical_bounds) = Xmag - X(end-1, horizontal_bounds, vertical_bounds);
elseif strcmp(wall,'west')
    X(1, horizontal_bounds, vertical_bounds) = Xmag - X(2, horizontal_bounds, vertical_bounds);
elseif strcmp(wall,'north')
    X( horizontal_bounds, end, vertical_bounds) = Xmag - X(horizontal_bounds, end-1, vertical_bounds);
elseif strcmp(wall,'south')
    X( horizontal_bounds, 1, vertical_bounds) = Xmag - X(horizontal_bounds, 2 ,vertical_bounds);
elseif strcmp(wall,'roof')
    X( horizontal_bounds, vertical_bounds, end) = Xmag - X(horizontal_bounds ,vertical_bounds, end-1);
elseif strcmp(wall,'floor')
    X( horizontal_bounds, vertical_bounds, 1) = Xmag - X(horizontal_bounds ,vertical_bounds, 2);
end

end