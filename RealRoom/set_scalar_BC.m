function Xi = set_scalar_BC(Xi,type)
%Homogeneous Dirichlet condition set on all walls of 3d box

if strcmp(type,'DDD')
Xieast = -Xi(end,:,:);
Xiwest = -Xi(1,:,:);
Xi = cat(1, Xiwest, Xi);
Xi = cat(1, Xi,Xieast);

Xinorth = -Xi(:,end,:);
Xisouth = -Xi(:,1,:);
Xi = cat(2,Xisouth, Xi);
Xi = cat(2, Xi, Xinorth);

Xifront = Xi(:,:,end);
Xiback =  Xi(:,:,1);
Xi = cat(3,Xiback,Xi);
Xi = cat(3, Xi, Xifront);
end
%Homogeneous Neumann condition for all walls in box'

if strcmp(type,'NNN')
Xieast = Xi(end,:,:);
Xiwest = Xi(1,:,:);
Xi = cat(1, Xiwest, Xi);
Xi = cat(1, Xi,Xieast);

Xinorth = Xi(:,end,:);
Xisouth = Xi(:,1,:);
Xi = cat(2,Xisouth, Xi);
Xi = cat(2, Xi, Xinorth);

Xifront = Xi(:,:,end);
Xiback =  Xi(:,:,1);
Xi = cat(3,Xiback,Xi);
Xi = cat(3, Xi, Xifront);
end
end

