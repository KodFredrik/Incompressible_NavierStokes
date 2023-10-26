function source = Add_2dflow_source(horizontal_bounds, vertical_bounds, strength)
%Creates a smooth velocity profile for inlets and outlets
%Outputs a 2-d array where values fall off smoothly from center in order to
%avoid sharp discontinuities between walls and inlets

lenH = size(horizontal_bounds);
lenV = size(vertical_bounds);
source = zeros(round(lenH(2)),round(lenV(2)));
H_midpoint = horizontal_bounds(round(lenH(2)/2));
V_midpoint = vertical_bounds(round(lenV(2)/2));
for i=horizontal_bounds
    for j = vertical_bounds
       if i == H_midpoint && j == V_midpoint
       source(i+1-horizontal_bounds(1),j+1-vertical_bounds(1)) = strength;
       else
       source(i+1-horizontal_bounds(1),j+1-vertical_bounds(1))  = strength / ...
           ((i-H_midpoint)^2 +(j-V_midpoint)^2);
       end
    end
end
end

