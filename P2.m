% 2D Chern insulator model, periodical boundary condition for both
% x-direction and y-direction.

Lx = 100;
Ly = 100;
m = 1;
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

kx = -pi+2*pi/Lx:2*pi/Lx:pi;
ky = -pi+2*pi/Ly:2*pi/Ly:pi;
% k = 0:2*pi/Lx:2*pi-2*pi/Lx;

for kxi = 1:length(kx)
    for kyi = 1:length(ky)
        H = sin(kx(kxi)).*sigma_x + sin(ky(kyi)).*sigma_y + (m+cos(kx(kxi))+cos(ky(kyi))).*sigma_z;
    
        e = eig(H);
        x = kx(kxi).*(zeros(2,1)+1);
        y = ky(kyi).*(zeros(2,1)+1);
        plot3(x,y,e,'.','color','k')
        hold on
    end
end