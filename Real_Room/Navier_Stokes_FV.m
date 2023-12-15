clear
close all
%Reynold's number and time step
V0 = 0.005;
L = 0.5;
Tmax = 20.1;
Tmin = 20;
Ximax = 0.000705; %Max measured value
g = 9.82;
alpha = 3.43 * 10^-3;
rho = 1.204 * 10^-6;
beta = 19 * 10^-6;
nu = 15.06 * 10^-6;
gamma = 16 * 10^-6;

Pr = nu/beta;
Re = V0*L/nu;
Fr = V0 / sqrt(g*L);
Sc = nu/gamma;
Ra = Pr*(g*alpha*(Tmax-Tmin)*L^3/nu^2);

Xi_const = 1/(Re*Sc);
T_const = 1/(Re*Pr);

radiator_real = 20.5;
window_real = 20;
window_pos = [15, 25];
radiator_pos = [15, 11];
window_size = 10;

%Grid Size N
Nx = 31;
Ny = 31;
Nz = 31;
dx_real = L/Nx;
dy_real = L/Ny;
dz_real = L/Nz;

dt_real = 0.001;
%Not including the boundaries, velocities are inbetween pressure points
Xi_real = zeros(Nx,Ny,Nz);
T_real = 0.5*(Tmax+Tmin)*ones(Nx,Ny,Nz);
U_real = zeros(Nx-1,Ny,Nz);
V_real = zeros(Nx,Ny-1,Nz);
W_real = zeros(Nx,Ny,Nz-1);

Q = zeros(Nx,Ny,Nz);
emission_rate = 0.000167; %Amount of CO2/s from one person
position = [16,16,16];
Q_real = Apply_source(Q, emission_rate, position, dt_real);

[dt,dx,dy,dz,U,V,W,T,Xi, Q, radiatorTemp, windowTemp] = Dimensionless_vars ...
    (dt_real,dx_real,dy_real,dz_real,U_real,V_real,W_real,T_real,Xi_real,...
    Q_real, radiator_real, window_real, V0,L,Tmin,Tmax, Ximax);

%Create linear initial temperature profile to speed up convergence
% for i=1:Nz
%     T(:,:,i) = FloorTemp - FloorTemp*i/Nz;
% end

%From floor to roof
% for i=1:Nz
%     T(:,:,i) = i/Nz * RoofTemp;
% end

%Number of iterations
timesteps = 1000000;
H = 1;
%Initialize plot
%[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
[Uplot,Vplot,Wplot] = stagger_back(U,V,W);
[f, h] = QuiverPlot(Uplot,Vplot,Wplot,Nx,Ny,Nz, H);

timediv = 1000;
movieVector1(timesteps/timediv) = struct('cdata',[],'colormap',[]);
%Outlet and inlet variables
Size = 10;
in_wall = 'west';
in_pos = [15,10];
out_wall = 'east';
out_pos = [15,25];

Umag_real = 0.005;
Umag = Umag_real / V0;
Vmag = 0;
Wmag = 0;

%Time iteration loop
time = dt;
for i = 1:timesteps

%Set Boundary Conditions
[Ubc,Vbc,Wbc] = Set_Room_BC(U,V,W);
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,Umag,Vmag,Wmag,in_pos,Size,in_wall);%inlet
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,-Umag,Vmag,Wmag,out_pos,Size,out_wall);%outlet
%Non-linear terms
[nonlinU,nonlinV,nonlinW] = nonlinear(Ubc,Vbc,Wbc,dx,dy,dz);
%Viscous terms
[viscousU, viscousV,viscousW] = viscous(Ubc,Vbc,Wbc, Re,dx,dy,dz);

Tavgz = average(T,3);
F_buoyancy = (Ra/(Re^2*Pr))* Tavgz;

%Compute intermediate velocities with explicit Euler, add force term here
Ustar = U + nonlinU.*dt + viscousU.*dt;
Vstar = V + nonlinV.*dt + viscousV.*dt;
Wstar = W + nonlinW.*dt + viscousW.*dt + F_buoyancy *dt - (1/Fr^2) *dt ;%

% %Solve Poisson's equation for pressure
[Ustarbc,Vstarbc,Wstarbc] = Set_Room_BC(Ustar,Vstar,Wstar);
[Ustarbc,Vstarbc,Wstarbc] = Add_Dir_part(Ustarbc,Vstarbc,Wstarbc,Umag,Vmag,Wmag,in_pos,Size,in_wall);%inlet
[Ustarbc,Vstarbc,Wstarbc] = Add_Dir_part(Ustarbc,Vstarbc,Wstarbc,-Umag,Vmag,Wmag,out_pos,Size,out_wall);%outlet

P = Solve_Poisson(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz,Nx,Ny,Nz,dt);

%Compute Unew, Vnew, Wnew
Px = diff(P,1,1)/dx;
Py = diff(P,1,2)/dy;
Pz = diff(P,1,3)/dz;

Uold = U;
Vold = V;
Wold = W;

U = Ustar - dt*Px;
V = Vstar - dt*Py;
W = Wstar - dt*Pz;

%Calculate Spread of scalar
[Ubc,Vbc,Wbc] = Set_Room_BC(U,V,W);
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,Umag,Vmag,Wmag,in_pos,Size,in_wall);%inlet
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,-Umag,Vmag,Wmag,out_pos,Size,out_wall);%outlet

Xi = Concentration_propagation(Xi,Xi_const, Ubc,Vbc,Wbc, Q,dt,dx,dy,dz,out_pos,Size,out_wall, in_pos, Size, in_wall);
T = Temperature_field(T, Ubc,Vbc,Wbc, T_const,windowTemp, radiatorTemp,window_pos, radiator_pos, window_size, dx,dy,dz,dt);

%Calculate some residuals to check for convergence
Uresidual = sqrt((U - Uold).^2);
Vresidual = sqrt((V - Vold).^2);
Wresidual = sqrt((W - Wold).^2);

TotResidual = sum(Uresidual, "all")/(Nz*Ny*(Nx-1)) + sum(Vresidual, "all")/(Nx*Nz*(Ny-1)) + sum(Wresidual, "all")/(Nx*Ny*(Nz-1));

disp("RMS Residual:")
disp(TotResidual)

%inlet and outlet regions
xrange = 1:10;
yrange = round(in_pos(1)-Size/2):round(in_pos(1)+Size/2);
zrange =  round(in_pos(2)-Size/2):round(in_pos(2)+Size/2);

Nr_elements = length(xrange)*length(yrange)*length(zrange);
UinletRes = Uresidual(xrange, yrange ,zrange);
VinletRes = Vresidual(xrange, yrange ,zrange);
WinletRes = Wresidual(xrange, yrange ,zrange);
InletResidual = sum(UinletRes, "all")/Nr_elements + sum(VinletRes, "all")/Nr_elements...
    + sum(WinletRes, "all")/Nr_elements;

disp("RMS inlet Residual:")
disp(InletResidual)

xrange = Nx-11:Nx-1;
yrange = round(out_pos(1)-Size/2):round(out_pos(1)+Size/2);
zrange =  round(out_pos(2)-Size/2):round(out_pos(2)+Size/2);

Nr_elements = length(xrange)*length(yrange)*length(zrange);
UoutletRes = Uresidual(xrange, yrange ,zrange);
VoutletRes = Vresidual(xrange, yrange ,zrange);
WoutletRes = Wresidual(xrange, yrange ,zrange);
OutletResidual = sum(UoutletRes, "all")/Nr_elements + sum(VoutletRes, "all")/Nr_elements...
    + sum(WoutletRes, "all")/Nr_elements;

disp("RMS outlet Residual:")
disp(OutletResidual)

if mod(i,timediv) == 0
    figure(1);
    title('Flow field at time = ', time*L/V0)
    [Uplot,Vplot,Wplot] = stagger_back(U,V,W);
    Uplot = permute(Uplot.*V0, [2,1,3]);
    Vplot = permute(Vplot.*V0, [2,1,3]);
    Wplot = permute(Wplot.*V0, [2,1,3]);

    h.UData = Uplot(1:H:Nx, 1:H:Ny,1:H:Nz);
    h.VData = Vplot(1:H:Nx, 1:H:Ny,1:H:Nz);
    h.WData = Wplot(1:H:Nx, 1:H:Ny,1:H:Nz);
    ColorQuiver(h);
    view([-60, 20])
    movieVector1(i/timediv) = getframe(f, [10 50 1500 1120]);
    pause(0.05)

end
    
if mod(i,timediv) == 0
    StreamSliceMovie(Uplot*V0,Vplot*V0,Wplot*V0,Nx,Ny,Nz,time*L/V0);
end

if mod(i,timediv) == 0
    ScalarMovie(Xi*Ximax,Nx,Ny,Nz,time*L/V0, Re, Sc, emission_rate)
end

if mod(i,timediv) == 0
    TemperatureMovie(T.*(Tmax-Tmin)+Tmin,Nx,Ny,Nz,time*L/V0, Pr, Re, Ra)
end

disp("Iteration:")
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
PressureMovie(P,Nx,Ny,Nz,time*L/V0)
