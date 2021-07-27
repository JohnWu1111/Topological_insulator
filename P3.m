% 2D Chern insulator model, periodical boundary condition for x-direction but
% finite size for y-direction, square lattice, after Fourier transform for
% x-axis.

tic;
Lx = 200;
Ly = 200;
m = -1;
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

% J = 1;

k = -pi+2*pi/Lx:2*pi/Lx:pi;
% k = 0:2*pi/Lx:2*pi-2*pi/Lx;

for ki = 1:length(k)
    H1 = zeros(Ly,Ly);
    H2 = zeros(Ly,Ly);
    H3 = zeros(Ly,Ly);
    for i = 1:Ly-1
        H1(i,i) = sin(k(ki));
        H3(i,i) = m + cos(k(ki));
        
        H2(i,i+1) = -1i/2;
        H2(i+1,i) = 1i/2;
        H3(i,i+1) = 1/2;
        H3(i+1,i) = 1/2;
    end
    H1(Ly,Ly) = sin(k(ki));
    H3(Ly,Ly) = m + cos(k(ki));
    
    H = kron(H1,sigma_x) + kron(H2,sigma_y) + kron(H3,sigma_z);
    
    e = eig(H);
    x = k(ki).*(zeros(Ly*2,1)+1);
    plot(x,e,'.','color','k')
    hold on
end
xlabel('kx')
ylabel('\epsilon_k')
str = strcat('m = ', num2str(m));
title(str)

toc;