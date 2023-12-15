function movieVectorC = concentrationMovie(Xi,Nx,Ny,Nz,time, Re, Sc, strength, movieVector2)
[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
pressureWriter = VideoWriter(['C:\Users\fredd\OneDrive\Skrivbord\flowfilms\Con_scatter_' num2str(time,'%.2f')], 'MPEG-4');
pressureWriter.FrameRate = 1;
%movieVector2(Nx) = struct('cdata',[],'colormap',[]);
%initialize pressureplot
f6 = figure(6);
f6.Position = [10 50 1040 800];
view(3)
%for k=1:Nx
    axis([1 Nx 1 Ny 1 Nz])
    perm = permute(Xi,[2,1,3]);
    x = [X(:); X(:); X(:)];
    y = [Y(:); Y(:); Z(:)];
    z = [Z(:); Z(:); Z(:)];
    s = [perm(:); perm(:); perm(:)];
    a
    S = scatter3(x,y,z, max(s, 1.e-15));
    axh = S(1).Parent;
    axis([1 Nx 1 Ny 1 Nz])
    title(axh, ['Concentration at time =' sprintf('%0.2f',time) ', Re = ' sprintf('%0.8f', Re) ...
        ', Sc = ' sprintf('%0.8f', Sc) ', Emission Rate = ' sprintf('%0.6f',strength)]) 
    xlabel(axh, 'X')
    ylabel(axh, 'Y')
    zlabel(axh, 'Z')
    movieVector2(k) = getframe(f6, [10 50 990 740]);
%end
open(pressureWriter);
writeVideo(pressureWriter,movieVector2);
close(pressureWriter);
end

