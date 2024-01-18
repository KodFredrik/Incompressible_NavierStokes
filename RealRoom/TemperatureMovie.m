function  TemperatureMovie(T,Nx,Ny,Nz,time, Pr, Re, Ra)
[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
pressureWriter = VideoWriter(['C:\Users\fredd\OneDrive\Skrivbord\flowfilms\Temperature_' num2str(time,'%.2f')], 'MPEG-4');
pressureWriter.FrameRate = 1;
movieVector2(Nx) = struct('cdata',[],'colormap',[]);
%initialize pressureplot
f5 = figure(5);
f5.Position = [10 50 1040 800];
view(3)
for k=1:Nx
    axis([1 Nx 1 Ny 1 Nz])
    S = slice(X,Y,Z,permute(T,[2,1,3]),[],k,[]);
    axh = S(1).Parent;
    axis([1 Nx 1 Ny 1 Nz])
    colorbar
    colormap jet
    title(axh, ['Temperature at time =' sprintf('%0.2f',time) ', Pr = ' sprintf('%0.8f', Pr) ', Re = ' sprintf('%0.5f',Re) ', Ra = ' sprintf('%0.5f',Ra)]) 
    xlabel(axh, 'X')
    ylabel(axh, 'Y')
    zlabel(axh, 'Z')
    movieVector2(k) = getframe(f5, [10 50 990 740]);
end
open(pressureWriter);
writeVideo(pressureWriter,movieVector2);
close(pressureWriter);
end

