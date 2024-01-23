function PressureMovie(P,Nx,Ny,Nz,time,path)
[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
pressureWriter = VideoWriter([path '\Pressure_' num2str(time,'%.2f')], 'MPEG-4');
pressureWriter.FrameRate = 1;
movieVector2(Nx) = struct('cdata',[],'colormap',[]);
%initialize pressureplot
f3 = figure(3);
f3.Position = [10 50 1040 800];
view(3)
for k=1:Nx
    S = slice(X,Y,Z,P,k,[],k);
    axh = S(1).Parent;
    colorbar
    colormap jet
    title(axh, 'Pressure at time =', time) 
    xlabel(axh, 'X')
    ylabel(axh, 'Y')
    zlabel(axh, 'Z')
    movieVector2(k) = getframe(f3, [10 50 990 740]);
end
open(pressureWriter);
writeVideo(pressureWriter,movieVector2);
close(pressureWriter);
end

