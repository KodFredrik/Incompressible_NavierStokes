function [XiUx, XiVy, XiWz] = Scalar_mixed_terms(U,V,W,Xi,dx,dy,dz)
%Compute mixed terms in concentration equation
XiAvgx = average(Xi,1);
XiAvgy = average(Xi,2);
XiAvgz = average(Xi,3);
UXi = U.*XiAvgx;
VXi = V.*XiAvgy;
WXi = W.*XiAvgz;

XiUx = diff(UXi,1,1)./dx;
XiVy = diff(VXi,1,2)./dy;
XiWz = diff(WXi,1,3)./dz;
end

