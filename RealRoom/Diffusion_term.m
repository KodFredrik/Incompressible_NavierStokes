function DiffXi = Diffusion_term(Xi, constant,dx,dy,dz)
%Compute diffusion term
Xixx= diff(Xi(:,2:end-1, 2:end-1), 2, 1);
Xiyy = diff(Xi(2:end-1,:, 2:end-1), 2, 2);
Xizz = diff(Xi(2:end-1, 2:end-1, :),2, 3);
DiffXi = constant.* (Xixx./(dx^2) + Xiyy./(dy^2) + Xizz./(dz^2));
end

