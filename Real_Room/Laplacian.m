function L = Laplacian(nx,ny,nz,dx,dy,dz)
%The function creates a three dimensional Laplacian with homogeneous
%Neumann conditions by using the kronecker product

%Begin with x-direction
Ix = ones(nx,1);
Dxx = spdiags([Ix -2*Ix Ix], [-1 0 1], nx, nx); %Creates 1-D Laplacian
%Sets Neumann condition for 1-D Laplacian
Dxx(1,1) = -1;
Dxx(1,2) = 1;
Dxx(end,end-1) = 1;
Dxx(end,end) = -1;
Dxx = Dxx / (dx^2);

%Same thing for y-direction
Iy = ones(ny,1);
Dyy = spdiags([Iy, -2*Iy Iy], [-1 0 1], ny, ny);
Dyy(1,1) = -1;
Dyy(1,2) = 1;
Dyy(end,end-1) = 1;
Dyy(end,end) = -1;
Dyy = Dyy/(dy^2);

%Same thing for z-direction
Iz = ones(nz,1);
Dzz = spdiags([Iz, -2*Iz Iz], [-1 0 1], nz, nz);
Dzz(1,1) = -1;
Dzz(1,2) = 1;
Dzz(end,end-1) = 1;
Dzz(end,end) = -1;
Dzz = Dzz / (dz^2);

%Computes kronecker product between 1-D laplacian and identity matrices to
%obtain 3-D Laplacian as a large sparse matrix.
L = kron(Dxx, kron(speye(ny),speye(nz)) ) + kron(speye(nx),kron(Dyy, speye(nz)) ) ...
    + kron(kron(speye(nx), speye(ny)), Dzz);

end

