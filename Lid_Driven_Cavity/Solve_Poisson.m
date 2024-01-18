function P = Solve_Poisson(Ustar,Vstar,Wstar,dx,dy,dz,Nx,Ny,Nz,dt)
%Solve Poisson's equation for pressure
[Ustarbc,Vstarbc,Wstarbc] = set_Dirichlet_BC(Ustar,Vstar,Wstar);
RHS = (1/dt) * div(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz);
L = Laplacian(Nx,Ny,Nz,dx,dy,dz);
%L(1,:) = 0;
%L(1,1) = 1;
%RHS = permute(RHS, [3,2,1]);
rhs = reshape(RHS,Nx*Ny*Nz, 1);
%rhs(1) = 0;

Pcorr = L\rhs;
P = reshape(Pcorr, [Nx,Ny,Nz]);
%P = permute(P, [3,2,1]);
end

