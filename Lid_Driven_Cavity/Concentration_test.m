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

load('C:\Users\fredd\OneDrive\Dokument\MATLAB\Lid_Driven_Cavity\Velocity_vars.mat')

Xi = zeros(Nx,Ny,Nz);
%Source term Q
Q = zeros(Nx,Ny,Nz);
Qstrength = 10000*0.000421;
position = [16,16,16];
Q = Apply_source(Q, Qstrength, position);
Diffusion_constant= 0.000016; %Co2 in air


%Number of iterations
timesteps = 100000;

%Time iteration loop
time = dt;
for i = 1:timesteps
Xi = Concentration_propagation(Xi,Diffusion_constant, U,V,W, Q,dt,dx,dy,dz);
if mod(i,10) == 0
    ScalarMovie(Xi,Nx,Ny,Nz,time, Diffusion_constant, Qstrength)
end

disp(i)
time = time + dt;
end