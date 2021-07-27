% Calculation of Chern number of a Chern insulator.
tic;
Lx = 1e4;
Ly = 1e4;
m = 1;

dkx = 2*pi/Lx;
dky = 2*pi/Ly;
kx = -pi:dkx:pi;
ky = -pi:dky:pi;

C = 0;
for i = 1:length(kx)-1
    for j = 1:length(ky)-1
        E = [sin(kx(i)) sin(ky(j)) m+cos(kx(i))+cos(ky(j))];
        n = E./norm(E);
        E_x = [sin(kx(i+1)) sin(ky(j)) m+cos(kx(i+1))+cos(ky(j))];
        n_x = E_x./norm(E_x);
        E_y = [sin(kx(i)) sin(ky(j+1)) m+cos(kx(i))+cos(ky(j+1))];
        n_y = E_y./norm(E_y);
        dnx = (n_x - n)./dkx;
        dny = (n_y - n)./dky;
        C = C + dot(n,cross(dnx,dny))*dkx*dky;
    end
end

C = C/(4*pi)
toc;