function [U,V,W] = stagger_back(U,V,W)

lenU = size(U);
Ueast = zeros(1,lenU(2),lenU(3));
Uwest = zeros(1,lenU(2),lenU(3));
U = cat(1, Uwest, U);
U = cat(1, U,Ueast);
U = average(U,1);

lenV = size(V);
Vnorth =  zeros(lenV(1),1,lenV(3));
Vsouth = zeros(lenV(1),1,lenV(3));
V = cat(2,Vsouth, V);
V = cat(2, V, Vnorth);
V = average(V,2);

lenW = size(W);
Wfront = zeros(lenW(1),lenW(2),1);
Wback =  zeros(lenW(1),lenW(2),1);
W = cat(3,Wback,W);
W = cat(3, W, Wfront);
W = average(W,3);
end

