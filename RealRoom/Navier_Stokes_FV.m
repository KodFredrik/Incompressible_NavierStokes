clear
close all

%Path to store video plots
path = 'C:\Users\fredd\OneDrive\Skrivbord\flowfilms';

%
g = 9.82;
alpha = 3.43 * 10^-3; %Volumetric expansion coefficient
rho = 1.204 * 10^-6; %Density of air
beta = 21.7 * 10^-6; %Thermal diffusivity of air
nu = 12.71 * 10^-6; %Kinematic viscosity of air
gamma = 16 * 10^-6; %Mass diffusion coefficient for CO2

%Characteristic quantities that determine the dimensionless quantities
V0 = 0.005;
L = 0.5;
Tmax = 20.1;
Tmin = 20;
Ximax = 0.000705; %Max measured value

%Dimensionless numbers
Pr = nu/beta;
Re = V0*L/nu;
Fr = V0 / sqrt(g*L);
Sc = nu/gamma;
Ra = Pr*(g*alpha*(Tmax-Tmin)*L^3/nu^2);

%Constants that go into the scalar propagation functions
Xi_const = 1/(Re*Sc);
T_const = 1/(Re*Pr);

%Number of cells N and step size in time and space
Nx = 31;
Ny = 31;
Nz = 31;
dx_real = L/Nx;
dy_real = L/Ny;
dz_real = L/Nz;
dt_real = 0.001;

%Hot and cold source position and temperature
%X,Y,Z,SIZE,TEMP
%Note that walls are given as 0 or N+1
temp_dir_array = [[0,20,20, 9, 20]; [Nx+1,10,10,9,20]; [16,0,16,13,20.1]; [10,Ny+1,10,11,20.1]];
temp_dir_array(:,5) = (temp_dir_array(:,5) - Tmin) ./ (Tmax-Tmin);

%Initialize arrays for all quantities
%The arrays do not include the boundaries, they will be concatenated to the
%array with a function.
%The different sizes for the velocities is because the grid is staggered,
%and they are defined between the pressure points.
Xi_real = zeros(Nx,Ny,Nz);
T_real = (0.3*Tmax+Tmin*0.7)*ones(Nx,Ny,Nz); %Temperature initialized with a temperature between the hot and cold source
U_real = zeros(Nx-1,Ny,Nz);
V_real = zeros(Nx,Ny-1,Nz);
W_real = zeros(Nx,Ny,Nz-1);
P = zeros(Nx,Ny,Nz);
Q = zeros(Nx,Ny,Nz);

%Define a source term that emits CO2 somewhere in the room
emission_rate = 0.000167; %Amount of CO2/s from one person
position1 = [16,16,16];
position2 = [10,10,10];
Q_1 = Apply_source(Q, emission_rate, position1);
Q_2 = Apply_source(Q, emission_rate, position2);
Q_real = Q_1 + Q_2;

%This function transforms all the real quantities to dimensionless form.
[dt,dx,dy,dz,U,V,W,T,Xi, Q] = Dimensionless_vars ...
    (dt_real,dx_real,dy_real,dz_real,U_real,V_real,W_real,T_real,Xi_real,...
    Q_real, V0,L,Tmin,Tmax, Ximax);


%Define how large the inlet and outlet should be, on which wall and where
%they should be and the air velocity at the inlet/outlet.
Size = 10; 
in_wall = 'west';
in_pos = [16,10]; 
out_wall = 'east';
out_pos = [16,25];

%This will be used for a dirichlet boundary condition.
Umag_real = 0.005;
Umag = Umag_real / V0; %Dimensionless transformation
Vmag = 0;
Wmag = 0;

%Number of iterations
timesteps = 1000;
%timediv determines how often a movie frame should be created
timediv = 500;

%Initialize vector field plot
H = 1; %Determines how many arrows should be in the vector field plot
[Uplot,Vplot,Wplot] = stagger_back(U,V,W);
[f, h] = QuiverPlot(Uplot,Vplot,Wplot,Nx,Ny,Nz, H);
movieVector1(timesteps/timediv) = struct('cdata',[],'colormap',[]);

%Initialize vectors to store values for Sensor 1 and set position
data_interval = timediv/10;
concentrations1 = zeros(timesteps/data_interval, 1);
sensortemps1 = zeros(timesteps/data_interval, 1);
sensor1x = 16;
sensor1y = 2;
sensor1z = 16;

%Same for sensor 2
concentrations2 = zeros(timesteps/data_interval, 1);
sensortemps2 = zeros(timesteps/data_interval, 1);
sensor2x = 26;
sensor2y = 30;
sensor2z = 26;


time = dt; %For plotting purposes
%The main loop
for i = 1:timesteps

%Set Boundary Conditions
[Ubc,Vbc,Wbc] = Set_Room_BC(U,V,W); %Sets no-slip on all walls
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,Umag,Vmag,Wmag,in_pos,Size,in_wall);%inlet
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,-Umag,Vmag,Wmag,out_pos,Size,out_wall);%standard outlet

%Experimenting with pressure condition
%[Ubc,Vbc,Wbc] = pressure_outlet_BC(Ubc,Vbc,Wbc,P,Vmag,Wmag,out_pos,Size,out_wall,dx,dy,dz, Re, dt);

%Compute non-linear terms
[nonlinU,nonlinV,nonlinW] = nonlinear(Ubc,Vbc,Wbc,dx,dy,dz);
%Compute viscous terms
[viscousU, viscousV,viscousW] = viscous(Ubc,Vbc,Wbc, Re,dx,dy,dz);

%Average temperature in z-direction to coincide with staggered W.
Tavgz = average(T,3);
%Compute Buoyant force with Boussinesq approximation
F_buoyancy = (Ra/(Re^2*Pr))* Tavgz;

%Compute intermediate velocities with explicit Euler, add force terms here
Ustar = U + nonlinU.*dt + viscousU.*dt;
Vstar = V + nonlinV.*dt + viscousV.*dt;
Wstar = W + nonlinW.*dt + viscousW.*dt + F_buoyancy *dt - (1/Fr^2) *dt ;

%The boundaries are added to the velocities once again because the
%divergence of the velocity field is computed in Poisson's equation
[Ustarbc,Vstarbc,Wstarbc] = Set_Room_BC(Ustar,Vstar,Wstar);
[Ustarbc,Vstarbc,Wstarbc] = Add_Dir_part(Ustarbc,Vstarbc,Wstarbc,Umag,Vmag,Wmag,in_pos,Size,in_wall);%inlet
[Ustarbc,Vstarbc,Wstarbc] = Add_Dir_part(Ustarbc,Vstarbc,Wstarbc,-Umag,Vmag,Wmag,out_pos,Size,out_wall);%standard outlet

%Experimenting
%[Ustarbc,Vstarbc,Wstarbc] = pressure_outlet_BC(Ustarbc,Vstarbc,Wstarbc,P,Vmag,Wmag,out_pos,Size,out_wall,dx,dy,dz, Re, dt);

%Solve Poisson's equation for pressure
P = Solve_Poisson(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz,Nx,Ny,Nz,dt);

%Compute Pressure gradient
Px = diff(P,1,1)/dx;
Py = diff(P,1,2)/dy;
Pz = diff(P,1,3)/dz;

%Store values at previous iteration in order to compute RMS-residual
Uold = U;
Vold = V;
Wold = W;

%Compute pressure correction to obtain new, divergence free velocity field.
U = Ustar - dt*Px;
V = Vstar - dt*Py;
W = Wstar - dt*Pz;

%Add boundary conditions to the velocity field required to solve the
%advection-diffusion equation
[Ubc,Vbc,Wbc] = Set_Room_BC(U,V,W);
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,Umag,Vmag,Wmag,in_pos,Size,in_wall);%inlet
[Ubc,Vbc,Wbc] = Add_Dir_part(Ubc,Vbc,Wbc,-Umag,Vmag,Wmag,out_pos,Size,out_wall);%standard outlet

%Experimenting
%[Ubc,Vbc,Wbc] = pressure_outlet_BC(Ubc,Vbc,Wbc,P,Vmag,Wmag,out_pos,Size,out_wall,dx,dy,dz, Re, dt);

%Calculate spreading of scalars (Concentration and temperature)
Xi = Concentration_propagation(Xi,Xi_const, Ubc,Vbc,Wbc, Q,dt,dx,dy,dz,out_pos,Size,out_wall, in_pos, Size, in_wall);
T = Temperature_field(T, Ubc,Vbc,Wbc, T_const, dx,dy,dz,dt, temp_dir_array);

%The main loop is done for the solution of the equations, the rest is just
%for calculation of residuals and plotting.

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

%Updates vector field
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

%Add frames to movie vectors
%Need to transform back from dimensionless quantities in input
if mod(i,timediv) == 0
    StreamSliceMovie(Uplot*V0,Vplot*V0,Wplot*V0,Nx,Ny,Nz,time*L/V0, path);
end

if mod(i,timediv) == 0
    ScalarMovie(Xi*Ximax,Nx,Ny,Nz,time*L/V0, Re, Sc, emission_rate, path)
end

if mod(i,timediv) == 0
    TemperatureMovie(T.*(Tmax-Tmin)+Tmin,Nx,Ny,Nz,time*L/V0, Pr, Re, Ra, path)
end

%Sensor1 plot data
if mod(i,data_interval) == 0
    concentrations1(i/data_interval) = Xi(sensor1x, sensor1y, sensor1z)*Ximax;
    sensortemps1(i/data_interval) = T(sensor1x, sensor1y,sensor1z).*(Tmax-Tmin)+Tmin;
end

%Sensor2 plot data
if mod(i,data_interval) == 0
    concentrations2(i/data_interval) = Xi(sensor2x, sensor2y, sensor2z)*Ximax;
    sensortemps2(i/data_interval) = T(sensor2x, sensor2y,sensor2z).*(Tmax-Tmin)+Tmin;
end

disp("Iteration:")
disp(i)
time = time + dt;
end
time = time - dt;

%Saving the movie
myWriter = VideoWriter([path '\flowfield'], 'MPEG-4');
myWriter.FrameRate = 1;
open(myWriter);
writeVideo(myWriter,movieVector1);
close(myWriter)

%Pressure at final time
PressureMovie(P,Nx,Ny,Nz,time*L/V0, path)

%Sensor1 plot
f10 = figure(10);
times = linspace(0,timesteps*dt_real, timesteps/data_interval);
hold on
plot(times, sensortemps1)
plot(times, sensortemps2)
hold off
title(['Sensor data at ', '[' sprintf('%0.1f', sensor1x), ',' sprintf('%0.1f', sensor1y), ',' sprintf('%0.1f', sensor1z) ']' ...
    ' and [' sprintf('%0.1f', sensor2x), ',' sprintf('%0.1f', sensor2y), ',' sprintf('%0.1f', sensor2z) ']'])
xlabel('Time [s]')
ylabel('Temperature')
legend('Sensor 1', 'Sensor 2')
saveas(f10, [path '\Sensor_temps.png'])

%Sensor2 plot
f11 = figure(11);
times = linspace(0,timesteps*dt_real, timesteps/data_interval);
hold on
plot(times, concentrations1)
plot(times, concentrations2)
hold off
title(['Sensor data at ', '[' sprintf('%0.1f', sensor1x), ',' sprintf('%0.1f', sensor1y), ',' sprintf('%0.1f', sensor1z) ']' ...
    ' and [' sprintf('%0.1f', sensor2x), ',' sprintf('%0.1f', sensor2y), ',' sprintf('%0.1f', sensor2z) ']'])
xlabel('Time [s]')
ylabel('Concentration')
legend('Sensor 1', 'Sensor 2')
saveas(f11, [path '\Sensor_concentrations.png'])
