function [U,V,W] = set_Dirichlet_BC(U,V,W)

%Dirichlet condition can be specified directly on east and west boundary
%for U. U on north, south, front or back is defined as 
%Unorth = ones(lenU(1),1,lenU(3))-U(:,end,:);

lenU = size(U);
Ueast = zeros(1,lenU(2),lenU(3));
Uwest = zeros(1,lenU(2),lenU(3));
U = cat(1, Uwest, U);
U = cat(1, U,Ueast);

Unorth = -U(:,end,:);
Usouth = -U(:,1,:);
U = cat(2,Usouth, U);
U = cat(2, U, Unorth);

lenU = size(U);
Ufront = 2*ones(lenU(1), lenU(2), 1)-U(:,:,end);
Uback =  -U(:,:,1);
U = cat(3,Uback,U);
U = cat(3, U, Ufront);

%Same thing for V

lenV = size(V);
Vnorth =  zeros(lenV(1),1,lenV(3));
Vsouth = zeros(lenV(1),1,lenV(3));
V = cat(2,Vsouth, V);
V = cat(2, V, Vnorth);

Veast = -V(end,:,:);
Vwest = -V(1,:,:);
V = cat(1, Vwest, V);
V = cat(1, V,Veast);

Vfront = -V(:,:,end);
Vback =  -V(:,:,1);
V = cat(3,Vback,V);
V = cat(3, V, Vfront);

%And W
lenW = size(W);
Wfront = zeros(lenW(1),lenW(2),1);
Wback =  zeros(lenW(1),lenW(2),1);
W = cat(3,Wback,W);
W = cat(3, W, Wfront);

Wnorth =  -W(:,end,:);
Wsouth = -W(:,1,:);
W = cat(2,Wsouth, W);
W = cat(2, W, Wnorth);

Weast = -W(end,:,:);
Wwest = -W(1,:,:);
W = cat(1, Wwest, W);
W = cat(1, W,Weast);

end

