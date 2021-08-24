% 2D Chern insulator model, periodical boundary condition for y-direction but
% finite size for x-direction, square lattice, after inversal Fourier transform for
% x-axis. Calculate the probability current of the edge state.

tic;
Lx = 100;
Ly = 100;
m = 1.5;
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

% J = 1;

k = -pi+2*pi/Ly:2*pi/Ly:pi;
% k = 0:2*pi/Lx:2*pi-2*pi/Lx;

phi_u = zeros(Ly,Ly);
phi_d = zeros(Ly,Ly);
    
phi_u_k = zeros(Ly,Ly);
phi_d_k = zeros(Ly,Ly);

phi_m_k = zeros(Lx,Ly);
phi_m = zeros(Lx,Ly);

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
    
    [phi,e] = eig(H);
    x = k(ki).*(zeros(Lx*2,1)+1);
    
    for j = 1:Lx       
%         phi_m_k(:,ki) = phi_m_k(:,ki) + phi(1:Lx,j);
        phi_m_k(:,ki) = phi_m_k(:,ki) + phi(Lx+1:end,j);
    end
    
%     plot(x,diag(e),'.','color','k')
%     hold on
end

for i = 1:Lx
    for j = 1:Ly
        for ki = 1:Ly
            phi_m(i,j) = phi_m(i,j) + phi_m_k(i,ki)*exp(1i*j*k(ki))/sqrt(Ly);
        end
    end
end

phi_m_c = conj(phi_m);
J_l = zeros(1,Ly);
J_r = zeros(1,Ly);
for i = 2:Ly-1
    J_l(i) = -1i*(phi_m_c(1,i)*(phi_m(1,i+1)-phi_m(1,i-1))/2 - phi_m(1,i)*(phi_m_c(1,i+1)-phi_m_c(1,i-1))/2);
    J_r(i) = -1i*(phi_m_c(end,i)*(phi_m(end,i+1)-phi_m(end,i-1))/2 - phi_m(end,i)*(phi_m_c(end,i+1)-phi_m_c(end,i-1))/2);
end
y = 2:Ly-1;
plot(y,J_l(2:Ly-1));
hold on
plot(y,J_r(2:Ly-1));

% xlabel('kx')
% ylabel('\epsilon_k')
% str = strcat('m = ', num2str(m));
% title(str)

toc;
