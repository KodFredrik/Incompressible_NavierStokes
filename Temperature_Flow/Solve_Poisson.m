function P = Solve_Poisson(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz,Nx,Ny,Nz,dt)
%Solve Poisson's equation for pressure
RHS = (1/dt) * div(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz);
L = Laplacian(Nx,Ny,Nz,dx,dy,dz);
%L(1,:) = 0;
%L(1,1) = 1;
rhs = reshape(RHS,Nx*Ny*Nz, 1);
%rhs(1) = 0;

Pcorr = L\rhs;
P = reshape(Pcorr, [Nx,Ny,Nz]);
end

