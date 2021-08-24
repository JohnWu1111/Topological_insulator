% 2D Chern insulator model, periodical boundary condition for y-direction but
% finite size for x-direction, square lattice, after inversal Fourier transform for
% x-axis.

tic;
Lx = 100;
Ly = 100;
m = 1.5;
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

% J = 1;

k = -pi+2*pi/Lx:2*pi/Lx:pi;
% k = 0:2*pi/Lx:2*pi-2*pi/Lx;

for ki = 1:length(k)
    H1 = zeros(Lx,Lx); % for sigma_x
    H2 = zeros(Lx,Lx); % for sigma_y
    H3 = zeros(Lx,Lx); % for sigma_z
    H4 = zeros(Lx,Lx); % for epsilon_k
    for i = 1:Lx-1
        H1(i,i) = sin(k(ki));
        H3(i,i) = m - cos(k(ki));
        H4(i,i) = -2*cos(k(ki));
        
        H2(i,i+1) = 1i/2;
        H2(i+1,i) = -1i/2;
        H3(i,i+1) = -1/2;
        H3(i+1,i) = -1/2;
        H4(i,i+1) = -1;
        H4(i+1,i) = -1;
    end
    H1(Lx,Lx) = sin(k(ki));
    H3(Lx,Lx) = m - cos(k(ki));
    H4(Lx,Lx) = -2*cos(k(ki));
    
    H = 3*(kron(H1,sigma_x) - kron(H2,sigma_y) + kron(H3,sigma_z)) + kron(H4,eye(2));
%     H = 3*(kron(H1,sigma_x) - kron(H2,sigma_y) + kron(H3,sigma_z));
    
    e = eig(H);
    x = k(ki).*(zeros(Lx*2,1)+1);
    plot(x,e,'.','color','k')
    hold on
end
xlabel('kx')
ylabel('\epsilon_k')
str = strcat('m = ', num2str(m));
title(str)

toc;
