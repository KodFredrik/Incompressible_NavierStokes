clear
close all
%Reynold's number and time step

dt = 0.001;
%Grid size N
Nx = 31;
Ny = 31;
Nz = 31;
dx = 1/Nx;
dy = 1/Ny;
dz = 1/Nz;

Re = 1000;

%Not including the boundaries, velocities are inbetween pressure points
Xi = zeros(Nx,Ny,Nz);
U = zeros(Nx-1,Ny,Nz);
V = zeros(Nx,Ny-1,Nz);
W = zeros(Nx,Ny,Nz-1);

%Source term Q
Q = zeros(Nx,Ny,Nz);
Qstrength = 10000*0.000421;
position = [15,15,15];
Q = Apply_source(Q, Qstrength, position);
Diffusion_constant= 0.000016; %Co2 in air

%Number of iterations
timesteps = 1000;
H = 3;
%Initialize plot
%[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
[Uplot,Vplot,Wplot] = stagger_back(U,V,W);
Uplot = permute(Uplot, [2,1,3]);
Vplot = permute(Vplot, [2,1,3]);
Wplot = permute(Wplot, [2,1,3]);
[f, h] = QuiverPlot(Uplot,Vplot,Wplot,Nx,Ny,Nz, H);

movieVector1(timesteps/1000) = struct('cdata',[],'colormap',[]);

%Time iteration loop
time = dt;
for i = 1:timesteps

%Set Boundary Conditions
[Ubc,Vbc,Wbc] = set_Dirichlet_BC(U,V,W);
%Non-linear terms
[nonlinU,nonlinV,nonlinW] = nonlinear(Ubc,Vbc,Wbc,dx,dy,dz);
%Viscous terms
[viscousU, viscousV,viscousW] = viscous(Ubc,Vbc,Wbc, Re,dx,dy,dz);

%Compute intermediate velocities with explicit Euler, add force term here
Ustar = U + nonlinU.*dt + viscousU.*dt;
Vstar = V + nonlinV.*dt + viscousV.*dt;
Wstar = W + nonlinW.*dt + viscousW.*dt;% - (g/Fr^2) * dt;

% %Solve Poisson's equation for pressure
P = Solve_Poisson(Ustar,Vstar,Wstar,dx,dy,dz,Nx,Ny,Nz,dt);

%Compute Unew, Vnew, Wnew
Px = diff(P,1,1)/dx;
Py = diff(P,1,2)/dy;
Pz = diff(P,1,3)/dz;

% divven = div(Ustar,Vstar,Wstar,dx,dy,dz);
% sum(divven, "all")

U = Ustar - dt*Px;
V = Vstar - dt*Py;
W = Wstar - dt*Pz;

% divven1 = div(U,V,W,dx,dy,dz);
% sum(divven1, "all")

%Calculate Spread of scalar
Xi = Concentration_propagation(Xi,Diffusion_constant, U,V,W, Q,dt,dx,dy,dz);
%

if mod(i,1000) == 0
    figure(1);
    title('Flow field at time = ', time)
    [Uplot,Vplot,Wplot] = stagger_back(U,V,W);
    Uplot = permute(Uplot, [2,1,3]);
    Vplot = permute(Vplot, [2,1,3]);
    Wplot = permute(Wplot, [2,1,3]);
    h.UData = Uplot(1:H:Nx, 1:H:Ny,1:H:Nz);
    h.VData = Vplot(1:H:Nx, 1:H:Ny,1:H:Nz);
    h.WData = Wplot(1:H:Nx, 1:H:Ny,1:H:Nz);
    ColorQuiver(h);
    view([-60, 20])
    movieVector1(i/1000) = getframe(f, [10 50 1500 1120]);
    pause(0.05)

end
    
if mod(i,1000) == 0
    StreamSliceMovie(Uplot,Vplot,Wplot,Nx,Ny,Nz,time);
end

if mod(i,1000) == 0
    ScalarMovie(Xi,Nx,Ny,Nz,time, Diffusion_constant, Qstrength)
end

disp(i)
time = time + dt;
end
time = time - dt;

%Saving the movie
myWriter = VideoWriter('C:\Users\fredd\OneDrive\Skrivbord\flowfilms\flowfield', 'MPEG-4');
myWriter.FrameRate = 1;
open(myWriter);
writeVideo(myWriter,movieVector1);
close(myWriter)

%Pressure at final time
PressureMovie(P,Nx,Ny,Nz,time)
