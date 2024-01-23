function Q = Apply_source(Q, emission_rate, position, dt)
%Q is an array of the same size as Xi, but contains only a distribution of
%an amount of concentration corresponding to the amount breathed out by a
%person during an interval dt
len = size(Q);
Q = zeros(len);
for i=1:len(1)
    for j=1:len(2)
        for k=1:len(3)
            if position == [i,j,k]
                Q(i,j,k) = 1;%emission_rate * dt;
            elseif abs(i-position(1))>5 || abs(j-position(2))>5 || abs(k-position(3)) > 5
                Q(i,j,k) = 0;
            else 
            Q(i,j,k) = 1 / ((i-position(1))^2 +(j-position(2))^2 + (k-position(3))^2);
            end
        end
    end
end
Q = emission_rate * Q ./ sum(Q, "all");
end

