function L = Laplacian(nx,ny,nz,dx,dy,dz)

Ix = ones(nx,1);
Dxx = spdiags([Ix -2*Ix Ix], [-1 0 1], nx, nx);
%Neumann condition
Dxx(1,1) = -1;
Dxx(1,2) = 1;
Dxx(end,end-1) = 1;
Dxx(end,end) = -1;
% Dxx(1,2) = 2;
% Dxx(end,end-1) = 2;
Dxx = Dxx / (dx^2);
%
Iy = ones(ny,1);
Dyy = spdiags([Iy, -2*Iy Iy], [-1 0 1], ny, ny);
%Neumann condition
Dyy(1,1) = -1;
Dyy(1,2) = 1;
Dyy(end,end-1) = 1;
Dyy(end,end) = -1;
% Dyy(1,2) = 2;
% Dyy(end,end-1) = 2;
Dyy = Dyy/(dy^2);
%
Iz = ones(nz,1);
Dzz = spdiags([Iz, -2*Iz Iz], [-1 0 1], nz, nz);
%Neumann condition
Dzz(1,1) = -1;
Dzz(1,2) = 1;
Dzz(end,end-1) = 1;
Dzz(end,end) = -1;
% Dyy(1,2) = 2;
% Dyy(end,end-1) = 2;
Dzz = Dzz / (dz^2);
%
L = kron(Dxx, kron(speye(ny),speye(nz)) ) + kron(speye(nx),kron(Dyy, speye(nz)) ) ...
    + kron(kron(speye(nx), speye(ny)), Dzz);
%Modify to avoid singular matrix. Matrix is not singular though?????
% L(1,:) = 0;
% L(1,1) = -2;

end

