function X = Set_Dirichlet_part(X,Xmag,pos,Size,wall)

horizontal_bounds = pos(1)-floor(Size/2):pos(1)+floor(Size/2);
vertical_bounds = pos(2)-floor(Size/2):pos(2)+floor(Size/2);
% 
Xmag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Xmag);

shapeEW = size(X(1, horizontal_bounds, vertical_bounds));
shapeNS = size(X( horizontal_bounds, 1, vertical_bounds));
shapeRF = size(X( horizontal_bounds, vertical_bounds, 1));

if strcmp(wall,'east')
    X(end, horizontal_bounds, vertical_bounds) = reshape(Xmag, shapeEW) - X(end-1, horizontal_bounds, vertical_bounds);
elseif strcmp(wall,'west')
    X(1, horizontal_bounds, vertical_bounds) = reshape(Xmag, shapeEW) - X(2, horizontal_bounds, vertical_bounds);
elseif strcmp(wall,'north')
    X( horizontal_bounds, end, vertical_bounds) = reshape(Xmag, shapeNS) - X(horizontal_bounds, end-1, vertical_bounds);
elseif strcmp(wall,'south')
    X( horizontal_bounds, 1, vertical_bounds) = reshape(Xmag, shapeNS) - X(horizontal_bounds, 2 ,vertical_bounds);
elseif strcmp(wall,'roof')
    X( horizontal_bounds, vertical_bounds, end) = reshape(Xmag, shapeRF) - X(horizontal_bounds, end-1 ,vertical_bounds);
elseif strcmp(wall,'floor')
    X( horizontal_bounds, vertical_bounds, 1) = reshape(Xmag, shapeRF) - X(horizontal_bounds, 2 ,vertical_bounds);
end
end

