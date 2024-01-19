function divergence = div(Ustarbc,Vstarbc,Wstarbc,dx,dy,dz)
%The function computes the divergence of the velocity field

%We are now calculating the divergence at the pressure point, which means
%adding BCs and computing a forward derivative on the staggered grid is
%equivalent to a central difference at the pressure point
Ubcx = diff(Ustarbc,1,1)/dx;
Vbcy = diff(Vstarbc,1,2)/dy;
Wbcz = diff(Wstarbc,1,3)/dz;
%Returns interior points
divergence = Ubcx(:,2:end-1,2:end-1) + Vbcy(2:end-1,:,2:end-1) ...
        +Wbcz(2:end-1,2:end-1,:);
end

