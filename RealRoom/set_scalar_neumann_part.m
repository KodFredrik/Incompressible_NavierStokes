function X = set_scalar_neumann_part(X,Xmag,pos,Size,wall,dx,dy,dz)
%Size is width of inlet, should be an odd number
%wall is 'east', 'west','north','south','roof','floor'

horizontal_bounds = pos(1)-floor(Size/2):pos(1)+floor(Size/2);
vertical_bounds = pos(2)-floor(Size/2):pos(2)+floor(Size/2);
% 
Xmag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Xmag);

%flow source returns array of wrong size, need to reshape
shapeEW = size(X(1, horizontal_bounds, vertical_bounds));
shapeNS = size(X( horizontal_bounds, 1, vertical_bounds));
shapeRF = size(X( horizontal_bounds, vertical_bounds, 1));

if strcmp(wall,'east')
    X(end, horizontal_bounds, vertical_bounds) = X(end-1, horizontal_bounds, vertical_bounds) + reshape(Xmag, shapeEW)*dx;
        %X(end-1, horizontal_bounds, vertical_bounds)*dx;
elseif strcmp(wall,'west')
    X(1, horizontal_bounds, vertical_bounds) = X(2, horizontal_bounds, vertical_bounds) - reshape(Xmag, shapeEW)*dx;
elseif strcmp(wall,'north')
    X( horizontal_bounds, end, vertical_bounds) = X( horizontal_bounds, end-1, vertical_bounds) + reshape(Xmag, shapeNS) *dy ;
elseif strcmp(wall,'south')
    X( horizontal_bounds, 1, vertical_bounds) = X(horizontal_bounds, 2 ,vertical_bounds) - reshape(Xmag, shapeNS) *dy;
elseif strcmp(wall,'roof')
    X( horizontal_bounds, vertical_bounds, end) = X(horizontal_bounds,vertical_bounds, end-1 ) + reshape(Xmag, shapeRF)*dz;
elseif strcmp(wall,'floor')
    X( horizontal_bounds, vertical_bounds, 1) = X(horizontal_bounds ,vertical_bounds, 2) - reshape(Xmag, shapeRF) *dz;
end

end