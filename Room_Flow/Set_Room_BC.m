function [U,V,W] = Set_Room_BC(U,V,W)
%Zero Dirichlet for the walls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_Ue = 0;
C_Uw = 0;
C_Un = 0;
C_Us = 0;
C_Ub = 0;
C_Uf = 0;

C_Ve = 0;
C_Vw = 0;
C_Vn = 0;
C_Vs = 0;
C_Vb = 0;
C_Vf = 0;

C_We = 0;
C_Ww = 0;
C_Wn = 0;
C_Ws = 0;
C_Wb = 0;
C_Wf = 0;

%Sets Dirichlet everywhere
lenU = size(U);
Ueast = C_Ue*ones(1,lenU(2),lenU(3));
Uwest = C_Uw*ones(1,lenU(2),lenU(3));
U = cat(1, Uwest, U);
U = cat(1, U,Ueast);

Unorth = C_Un*2*ones(lenU(1)+2, 1, lenU(3)) - U(:,end,:);
Usouth = C_Us*2*ones(lenU(1)+2, 1, lenU(3)) - U(:,1,:);
U = cat(2,Usouth, U);
U = cat(2, U, Unorth);

Ufront = C_Uf*2*ones(lenU(1)+2, lenU(2)+2, 1)-U(:,:,end);
Uback =  C_Ub*2*ones(lenU(1)+2, lenU(2)+2, 1)-U(:,:,1);
U = cat(3,Uback,U);
U = cat(3, U, Ufront);

%Same thing for V
lenV = size(V);
Vnorth =  C_Vn*ones(lenV(1),1,lenV(3));
Vsouth = C_Vs*ones(lenV(1),1,lenV(3));
V = cat(2,Vsouth, V);
V = cat(2, V, Vnorth);

Veast = C_Ve*2*ones(1, lenV(2)+2, lenV(3)) - V(end,:,:);
Vwest = C_Vw*2*ones(1, lenV(2)+2, lenV(3)) - V(1,:,:);
V = cat(1, Vwest, V);
V = cat(1, V,Veast);

Vfront = C_Vf*2*ones(lenV(1)+2, lenV(2)+2, 1) - V(:,:,end);
Vback = C_Vb*2*ones(lenV(1)+2, lenV(2)+2, 1) - V(:,:,1);
V = cat(3,Vback,V);
V = cat(3, V, Vfront);

%And W
lenW = size(W);
Wfront = C_Wf*ones(lenW(1),lenW(2),1);
Wback =  C_Wb*ones(lenW(1),lenW(2),1);
W = cat(3,Wback,W);
W = cat(3, W, Wfront);

Wnorth = C_Wn*2*ones(lenW(1), 1 , lenW(3)+2) - W(:,end,:);
Wsouth = C_Ws*2*ones(lenW(1), 1 , lenW(3)+2)  - W(:,1,:);
W = cat(2,Wsouth, W);
W = cat(2, W, Wnorth);

Weast = C_We*2*ones(1, lenW(2)+2, lenW(3)+2) - W(end,:,:);
Wwest = C_Ww*2*ones(1, lenW(2)+2, lenW(3)+2) - W(1,:,:);
W = cat(1, Wwest, W);
W = cat(1, W,Weast);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

