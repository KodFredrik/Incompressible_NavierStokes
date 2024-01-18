function [U,V,W] = Add_neumann_part(U,V,W,Umag,Vmag,Wmag,pos,Size,wall,dx,dy,dz)
%U,V,W is expected to already have zero Dirichlet BC, i.e. walls.
%Umag is intensity of air velocity
%Upos is two coordinates to denote center of inlet placement
%Size is width of inlet, should be an odd number
%wall is 'east', 'west','north','south','roof','floor'

horizontal_bounds = pos(1)-floor(Size/2):pos(1)+floor(Size/2);
vertical_bounds = pos(2)-floor(Size/2):pos(2)+floor(Size/2);
% 
Umag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Umag);
Vmag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Vmag);
Wmag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Wmag);

%flow source returns array of wrong size, need to reshape
shapeEW = size(U(1, horizontal_bounds, vertical_bounds));
shapeNS = size(U( horizontal_bounds, 1, vertical_bounds));
shapeRF = size(U( horizontal_bounds, vertical_bounds, 1));

if strcmp(wall,'east')
    U(end, horizontal_bounds, vertical_bounds) = U(end-2, horizontal_bounds, vertical_bounds) + 2*reshape(Umag,shapeEW)*dx;
    V(end, horizontal_bounds, vertical_bounds) = V(end-1, horizontal_bounds, vertical_bounds) + reshape(Vmag, shapeEW)*dx;
    W(end, horizontal_bounds, vertical_bounds) =  W(end-1, horizontal_bounds, vertical_bounds) + reshape(Wmag,shapeEW)*dx;
elseif strcmp(wall,'west')
    U(1, horizontal_bounds, vertical_bounds) = U(3, horizontal_bounds, vertical_bounds) - 2*reshape(Umag,shapeEW)*dx;
    V(1, horizontal_bounds, vertical_bounds) = V(2, horizontal_bounds, vertical_bounds) - reshape(Vmag, shapeEW)*dx;
    W(1, horizontal_bounds, vertical_bounds) = W(2, horizontal_bounds, vertical_bounds) - reshape(Wmag, shapeEW)*dx ;
elseif strcmp(wall,'north')
    U( horizontal_bounds, end, vertical_bounds) = U( horizontal_bounds, end-1, vertical_bounds) + reshape(Umag, shapeNS) *dy ;
    V( horizontal_bounds, end, vertical_bounds) = V( horizontal_bounds, end-2, vertical_bounds) + 2* reshape(Vmag, shapeNS) *dy;
    W( horizontal_bounds, end, vertical_bounds) = W( horizontal_bounds, end-1, vertical_bounds) + reshape(Wmag,shapeNS) *dy;
elseif strcmp(wall,'south')
    U( horizontal_bounds, 1, vertical_bounds) = U(horizontal_bounds, 2 ,vertical_bounds) - reshape(Umag, shapeNS) *dy;
    V( horizontal_bounds, 1, vertical_bounds) = V( horizontal_bounds, 3, vertical_bounds) - 2* reshape(Vmag, shapeNS) *dy;
    W( horizontal_bounds, 1, vertical_bounds) = W( horizontal_bounds, 2, vertical_bounds) - reshape(Wmag, shapeNS) *dy;
elseif strcmp(wall,'roof')
    U( horizontal_bounds, vertical_bounds, end) = U(horizontal_bounds,vertical_bounds, end-1 ) + reshape(Umag, shapeRF)*dz;
    V( horizontal_bounds, vertical_bounds, end) = V(horizontal_bounds ,vertical_bounds, end-1) + reshape(Vmag, shapeRF) *dz;
    W( horizontal_bounds, vertical_bounds, end) = W( horizontal_bounds, vertical_bounds, end-2)+ 2* reshape(Wmag, shapeRF)*dz;
elseif strcmp(wall,'floor')
    U( horizontal_bounds, vertical_bounds, 1) = U(horizontal_bounds ,vertical_bounds, 2) - reshape(Umag, shapeRF) *dz;
    V( horizontal_bounds, vertical_bounds, 1) = V(horizontal_bounds ,vertical_bounds, 2) - reshape(Vmag, shapeRF) *dz;
    W( horizontal_bounds, vertical_bounds, 1) = W( horizontal_bounds, vertical_bounds, 3) - 2* reshape(Wmag, shapeRF) *dz;
end

end
