
%Script to investigate expression from von neumann analysis of
%advection-diffusion equation with central diffrence to first derivative
%and central difference to second derivative




f = figure(1);

%dx = linspace(0, 0.2, 100);
%dt = linspace(0, 0.1, 100);
dtstart = 0;
dtend = 0.1;
dtincr = dtend/50;

[ dt , x] = meshgrid( dtstart: dtincr : dtend, -pi:0.1:pi);
%Re = linspace(0, 10000, 100);
%Re = 1000;


u0 = 0.1;
nu = 16.82 * 10^-6;
L = 3;
dx = L/130;
Re = u0*L/nu;
%Re = 17830;
c = 1/Re;
%x = linspace(-pi, pi, 101);

%Gmod = 1 + (dt./dx).^2 .*sin(x).^2 - (4.*dt./(dx.^2.*Re) + 4.*dt.^2./(dx.^4.*Re.^2));
a = u0.*dt./dx;
%b = 4.*dt.*c./dx.^2; %this is for central difference on advective
b = (4*c.*dt./(dx.^2)) - (2*u0.*dt./dx);
Gmod = 1 + a.^2 .* sin(x).^2 - 2.*b.*sin(x/2).^2 + b.^2 .* sin(x/2).^4;
f = surf(dt, x, Gmod);


len = size(dt);
for t = 1:len(2)

    if max(Gmod(:,t)) > 1
        disp([num2str(dt(t-1, t-1)) "is largerst time step for stability"])
        title(["largest time step for stability: " + num2str(dt(t-1, t-1)) + "s" , "Re = " + num2str(Re) ...
            " u0 = " + num2str(u0) , "dx = " + num2str(dx)]);
        break
    end
end

xlabel('dt')
ylabel('phi')
zlabel('|G|^2')
disp("I hope you had a good loop")

