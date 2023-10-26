function [f, h] = QuiverPlot(Uplot,Vplot,Wplot,Nx,Ny,Nz, H)
%Initialize plot
[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
f = figure(1);
f.Position = [10 50 1560 1200];
%H = 2;
h = quiver3(X(1:H:Nx, 1:H:Ny,1:H:Nz),Y(1:H:Nx, 1:H:Ny,1:H:Nz),Z(1:H:Nx, 1:H:Ny,1:H:Nz)...
    ,Uplot(1:H:Nx, 1:H:Ny,1:H:Nz),Vplot(1:H:Nx, 1:H:Ny,1:H:Nz),Wplot(1:H:Nx, 1:H:Ny,1:H:Nz));
set(h,'autoscale','on','autoscalefactor',1.5);
grid on
colormap jet
colorbar
%
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
end

