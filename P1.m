% 2D lattice model, periodical boundary condition for x-direction but
% finite size for y-direction, square lattice, after Fourier transform for
% x-axis.

Lx = 100;
Ly = 100;
J = 1;

k = -pi+2*pi/Lx:2*pi/Lx:pi;
% k = 0:2*pi/Lx:2*pi-2*pi/Lx;

for ki = 1:length(k)
    H = zeros(Ly,Ly);
    for i = 1:Ly-1
        H(i,i+1) = -J;
        H(i+1,i) = -J;
        H(i,i) = -2*J*cos(k(ki));
    end
    H(Ly,Ly) = -2*J*cos(k(ki));
    
    e = eig(H)
    x = k(ki).*(zeros(Ly,1)+1);
    plot(x,e,'.','color','k')
    hold on
end