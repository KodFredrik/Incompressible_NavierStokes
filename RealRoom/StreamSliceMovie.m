function StreamSliceMovie(Uplot,Vplot,Wplot,Nx,Ny,Nz,time, path)
%Streamslicemovie
streamsliceWriter = VideoWriter([path '\Streamslice_'  num2str(time,'%.2f')], 'MPEG-4');
streamsliceWriter.FrameRate = 1;
movieVector3(Nx+Ny+Nz) = struct('cdata',[],'colormap',[]);
%Calculate magnitude of velocity
VelocityMagnitude = sqrt(Uplot.^2+Vplot.^2+Wplot.^2);
[X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz);
f2 = figure(2);
f2.Position = [10 50 1040 800];
for k=1:Nx
    clf
    view(3)
    [verts, averts] = streamslice(X,Y,Z,Uplot,Vplot,Wplot,k, [], []);
    %[verts, averts] = streamslice(X,Y,Z,fliplr(Uplot),fliplr(Vplot),fliplr(Wplot),k, [], []);
    slice(X,Y,Z,VelocityMagnitude,k,[],[]);
    streamline([verts, averts]);
    shading interp;
    axis([1 Nx 1 Ny 1 Nz])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Streamlines at time = ', time)
    colorbar
    colormap jet
    movieVector3(k) = getframe(f2, [10 50 990 740]);
end
for i=1:Ny
    clf
    view(3)
    [verts, averts] = streamslice(X,Y,Z,Uplot,Vplot,Wplot,[], i, []);
    %[verts, averts] = streamslice(X,Y,Z,fliplr(Uplot),fliplr(Vplot),fliplr(Wplot),[], i, []);
    slice(X,Y,Z,VelocityMagnitude,[],i,[]);
    streamline([verts, averts]);
    shading interp;
    axis([1 Nx 1 Ny 1 Nz])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Streamlines at time = ', time)
    colorbar
    colormap jet
    movieVector3(Nx+i) = getframe(f2, [10 50 990 740]);
end
for j=1:Nz
    clf
    view(3)
    [verts, averts] = streamslice(X,Y,Z,Uplot,Vplot,Wplot,[], [], j);
    %[verts, averts] = streamslice(X,Y,Z,fliplr(Uplot),fliplr(Vplot),fliplr(Wplot),[], [], j);
    slice(X,Y,Z,VelocityMagnitude,[],[],j);
    streamline([verts, averts]);
    shading interp;
    axis([1 Nx 1 Ny 1 Nz])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Streamlines at time = ', time)
    colorbar
    colormap jet
    movieVector3(Nx+Ny+j) = getframe(f2, [10 50 990 740]);
end
open(streamsliceWriter);
writeVideo(streamsliceWriter,movieVector3);
close(streamsliceWriter);
end

