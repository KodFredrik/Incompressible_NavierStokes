function Xi = Apply_source(Xi, strength, position)
len = size(Xi);

for i=1:len(1)
    for j=1:len(2)
        for k=1:len(3)
            if position == [i,j,k]
                Xi(i,j,k) = strength;
            elseif abs(i-position(1))>5 || abs(j-position(2))>5 || (k-position(3)) > 5
                Xi(i,j,k) = 0;
            else 
            Xi(i,j,k) = strength / ((i-position(1))^2 ... 
                +(j-position(2))^2 + (k-position(3))^2);
            end
        end
    end
end
end

