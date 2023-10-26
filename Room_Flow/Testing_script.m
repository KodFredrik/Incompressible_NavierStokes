% [a, b, c] = meshgrid(1:7,1:7,1:7);
% A = 0*ones(5,5,5);
% B = 0*ones(5,5,5);
% C = 0*ones(5,5,5);
% [A,B,C] = Set_Room_BC(A,B,C);
% Umag = 0;
% Vmag = 0;
% Wmag = 0;
% pos = [3,3];
% size = 3;
% wall = 'east';
% dx = 1;
% dy = 1;
% dz = 1;
% 
% [A,B,C] = Add_outlet(A,B,C,Umag,Vmag,Wmag,pos,size,wall,dx,dy,dz);
% %[A,B,C] = stagger_back(Abc,Bbc,Cbc);
% figure(3)
% abc = quiver3(a,b,c,permute(A,[2,1,3]),permute(B,[2,1,3]),permute(C,[2,1,3]));
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

horizontal_bounds = 6:14;
vertical_bounds = 6:14;
pos = [9,9];
strength = 2;

lenH = size(horizontal_bounds);
lenV = size(vertical_bounds);



source = Add_2dflow_source(horizontal_bounds, vertical_bounds,strength);
