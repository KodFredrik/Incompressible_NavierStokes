function P = Solve_Poisson(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz,Nx,Ny,Nz,dt)
%The function solves Poisson's equation for pressure and returns the
%pressure P.

%The right hand side is computed
RHS = (1/dt) * div(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz);
%The Laplacian is defined
L = Laplacian(Nx,Ny,Nz,dx,dy,dz);

rhs = reshape(RHS,Nx*Ny*Nz, 1);
Pcorr = L\rhs;
P = reshape(Pcorr, [Nx,Ny,Nz]);

end

