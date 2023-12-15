function P = Solve_Poisson(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz,Nx,Ny,Nz,dt)
%Solve Poisson's equation for pressure
RHS = (1/dt) * div(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz);
L = Laplacian(Nx,Ny,Nz,dx,dy,dz);

%RHS = permute(RHS, [3,2,1]);
rhs = reshape(RHS,Nx*Ny*Nz, 1);
Pcorr = L\rhs;


% perp = symamd(L);
% Rp = chol(L(perp,perp));
% Rpt = Rp';
% p(perp) = -Rp\(Rpt\RHS(perp));
% P = reshape(p, [Nx,Ny,Nz]);

P = reshape(Pcorr, [Nx,Ny,Nz]);
%P = permute(P, [3,2,1]);
end

