clear
close all
dt = 0.001;
%Grid size N
Nx = 31;
Ny = 31;
Nz = 31;
dx = 1/Nx;
dy = 1/Ny;
dz = 1/Nz;

load('C:\Users\fredd\OneDrive\Dokument\MATLAB\Room_flow\Velocity_vars,.mat')

Xi = zeros(Nx,Ny,Nz);
%Source term Q
Q = zeros(Nx,Ny,Nz);
Qstrength =10000* 0.000421;
position = [16,16,16];
Q = Apply_source(Q, Qstrength, position);
Diffusion_constant= 1000*0.000016; %Co2 in air

out_wall = 'east';
out_pos = [20,20];
Size = 9;
in_wall = 'west';
in_pos = [10,10];

%Number of iterations
timesteps = 100000;

%Time iteration loop
time = dt;
for i = 1:timesteps

[Ubc,Vbc,Wbc] = Set_Room_BC(U,V,W);
[Ubc,Vbc,Wbc] = Add_inlet(Ubc,Vbc,Wbc,2,0,0,in_pos,Size,in_wall);
[Ubc,Vbc,Wbc] = Add_outlet(Ubc,Vbc,Wbc,0,0,0,out_pos,Size,out_wall,dx,dy,dz);
Xi = Concentration_propagation(Xi,Diffusion_constant, Ubc,Vbc,Wbc, Q,dt,dx,dy,dz,out_pos,Size,out_wall);
if mod(i,1000) == 0
    ScalarMovie(Xi,Nx,Ny,Nz,time, Diffusion_constant, Qstrength)
end

disp(i)
time = time + dt;
end