function  ScalarMovie(Xi,Nx,Ny,Nz,time, Re, Sc, strength)
[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
pressureWriter = VideoWriter(['C:\Users\fredd\OneDrive\Skrivbord\flowfilms\Concentration_' num2str(time,'%.2f')], 'MPEG-4');
pressureWriter.FrameRate = 1;
movieVector2(Nx) = struct('cdata',[],'colormap',[]);
%initialize pressureplot
f4 = figure(4);
f4.Position = [10 50 1040 800];
view(3)
for k=1:Nx
    axis([1 Nx 1 Ny 1 Nz])
    S = slice(X,Y,Z,permute(Xi,[2,1,3]),[],k,[]);
    axh = S(1).Parent;
    axis([1 Nx 1 Ny 1 Nz])
    colorbar
    colormap jet
    title(axh, ['Concentration at time =' sprintf('%0.2f',time) ', Re = ' sprintf('%0.8f', Re) ...
        ', Sc = ' sprintf('%0.8f', Sc) ', Emission Rate = ' sprintf('%0.6f',strength)]) 
    xlabel(axh, 'X')
    ylabel(axh, 'Y')
    zlabel(axh, 'Z')
    movieVector2(k) = getframe(f4, [10 50 990 740]);
end
open(pressureWriter);
writeVideo(pressureWriter,movieVector2);
close(pressureWriter);
end
