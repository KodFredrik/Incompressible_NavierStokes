
function avgA = average(A, dim)
%Averages a 3d array along some dimension. Used to calculate cross terms.
    if(dim==1)
        avgA = 0.5*(A(2:end,:,:)+A(1:end-1,:,:));
    elseif(dim==2)
        avgA = 0.5*(A(:,2:end,:)+A(:,1:end-1,:));
    elseif(dim==3)
        avgA = 0.5*(A(:,:,2:end)+A(:,:,1:end-1));
    end
end