function [U,V,W] = Add_Dir_part(U,V,W,Umag,Vmag,Wmag,pos,Size,wall)
%This function sets a dirichlet condition at a square section of the wall.
%The function is intended to mimic the velocity profile from a pipe, and
%therefore creates a quadratic profile for the velocities.

%U,V,W is expected to already have zero Dirichlet BC, i.e. walls.
%Umag is intensity of air velocity
%Upos is two coordinates to denote center of inlet/outlet placement
%Size is width of inlet/outlet, preferably an odd number because because
%the center point of an even number is not an integer.
%wall is 'east', 'west','north','south','roof','floor'

%Rounding uneven numbers, for example, 4 and 5 will give inlet of same
%Size. 
horizontal_bounds = pos(1)-floor(Size/2):pos(1)+floor(Size/2);
vertical_bounds = pos(2)-floor(Size/2):pos(2)+floor(Size/2);

%This creates a quadratic velocity profile
Umag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Umag);
Vmag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Vmag);
Wmag = Add_2dflow_source(horizontal_bounds, vertical_bounds, Wmag);

%Get shape to be able to add velocity profile to array
shapeEW = size(U(1, horizontal_bounds, vertical_bounds));
shapeNS = size(U( horizontal_bounds, 1, vertical_bounds));
shapeRF = size(U( horizontal_bounds, vertical_bounds, 1));

%Sets dirichlet condition. Note that a negative value for Umag will give an
%outlet and a positive will give an inlet.
if strcmp(wall,'east')
    U(end, horizontal_bounds, vertical_bounds) = -reshape(Umag, shapeEW);
    V(end, horizontal_bounds, vertical_bounds) = reshape(Vmag, shapeEW) - V(end-1, horizontal_bounds, vertical_bounds);
    W(end, horizontal_bounds, vertical_bounds) = reshape(Wmag, shapeEW) - W(end-1, horizontal_bounds, vertical_bounds);
elseif strcmp(wall,'west')
    U(1, horizontal_bounds, vertical_bounds) = reshape(Umag, shapeEW);
    V(1, horizontal_bounds, vertical_bounds) = reshape(Vmag, shapeEW) - V(2, horizontal_bounds, vertical_bounds);
    W(1, horizontal_bounds, vertical_bounds) = reshape(Wmag, shapeEW) - W(2, horizontal_bounds, vertical_bounds);
elseif strcmp(wall,'north')
    U( horizontal_bounds, end, vertical_bounds) = reshape(Umag, shapeNS) - U( horizontal_bounds, end-1, vertical_bounds);
    V( horizontal_bounds, end, vertical_bounds) = - reshape(Vmag, shapeNS);
    W( horizontal_bounds, end, vertical_bounds) = reshape(Wmag, shapeNS) - W( horizontal_bounds, end-1, vertical_bounds);
elseif strcmp(wall,'south')
    U( horizontal_bounds, 1, vertical_bounds) = reshape(Umag, shapeNS) - U(horizontal_bounds, 2 ,vertical_bounds);
    V( horizontal_bounds, 1, vertical_bounds) = reshape(Vmag, shapeNS);
    W( horizontal_bounds, 1, vertical_bounds) = reshape(Wmag, shapeNS) - W( horizontal_bounds, 2, vertical_bounds);
elseif strcmp(wall,'roof')
    U( horizontal_bounds, vertical_bounds, end) = reshape(Umag, shapeRF) - U(horizontal_bounds,vertical_bounds, end-1 );
    V( horizontal_bounds, vertical_bounds, end) = reshape(Vmag, shapeRF) - V(horizontal_bounds ,vertical_bounds, end-1);
    W( horizontal_bounds, vertical_bounds, end) = -reshape(Wmag, shapeRF);
elseif strcmp(wall,'floor')
    U( horizontal_bounds, vertical_bounds, 1) = reshape(Umag, shapeRF) - U(horizontal_bounds,vertical_bounds, 2 );
    V( horizontal_bounds, vertical_bounds, 1) = reshape(Vmag, shapeRF) - V(horizontal_bounds,vertical_bounds, 2 );
    W( horizontal_bounds, vertical_bounds, 1) = reshape(Wmag, shapeRF);
end

end

