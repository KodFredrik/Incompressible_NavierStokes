
%function source = Add_2dflow_source(horizontal_bounds, vertical_bounds, strength)
%Creates a smooth velocity profile for inlets and outlets
%Outputs a 2-d array where values fall off smoothly from center in order to
%avoid sharp discontinuities between walls and inlets

strength = 5;
horizontal_bounds = 2:16;
vertical_bounds = 2:16;
[X, Y] = meshgrid(horizontal_bounds, vertical_bounds);

lenH = length(horizontal_bounds); %size can be replaced by length to clean up a bit...
lenV = length(vertical_bounds);
source = zeros(round(lenH),round(lenV));
H_midpoint = horizontal_bounds(round(lenH/2));
V_midpoint = vertical_bounds(round(lenV/2));

for i=horizontal_bounds
    for j = vertical_bounds
       dist = 1 - ((i-H_midpoint) / (0.5*lenH) )^2 - ( (j-V_midpoint) / (0.5*lenV) )^2; %distance from origin
       if i == H_midpoint && j == V_midpoint
       source(i+1-horizontal_bounds(1),j+1-vertical_bounds(1)) = strength;
       elseif dist < 0 
           source(i+1-horizontal_bounds(1),j+1-vertical_bounds(1)) = 0;
       else  
       %source(i+1-horizontal_bounds(1),j+1-vertical_bounds(1))  = strength
       %/ sqrt(((i-H_midpoint)^2 +(j-V_midpoint)^2)); % falls off as 1/r
       source(i+1-horizontal_bounds(1),j+1-vertical_bounds(1)) =  strength * dist;
       end
    end
end

surf(X,Y,source)
