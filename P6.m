% 2D Chern insulator model, periodical boundary condition for y-direction but
% finite size for x-direction, square lattice, after inversal Fourier transform for
% x-axis. Calculate the probability current of the edge state.

tic;
Lx = 40;
Ly = 40;
m = 1.5;
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

L = Lx*Ly;
H1 = zeros(L,L); % for sigma_x
H2 = zeros(L,L); % for sigma_y
H3 = zeros(L,L); % for sigma_z
H4 = zeros(L,L); % for epsilon_k

% diagonal term
for i = 1:Lx
    for j = 1:Ly
        k = j + (i-1)*Lx;
        H3(k,k) = m;
    end
end

% off-diagonal term, x
for j = 1:Ly
    for i = 1:Lx-1    
        k = j + (i-1)*Lx;
        H2(k,k+Lx) = -1i/2;
        H2(k+Lx,k) = 1i/2;
        H3(k,k+Lx) = -1/2;
        H3(k+Lx,k) = -1/2;
        H4(k,k+Lx) = -1;
        H4(k+Lx,k) = -1;
    end
%     k = Lx + (j-1)*Ly;
%     k0 = (j-1)*Ly;
%     H2(k,k0) = 1i/2;
%     H2(k0,k) = -1i/2;
%     H3(k,k0) = 1/2;
%     H3(k0,k) = 1/2;
%     H4(k,k0) = 1;
%     H4(k0,k) = 1;
end

% off-diagonal term, y
for i = 1:Lx
    for j = 1:Ly-1    
        k = j + (i-1)*Lx;
        H1(k,k+1) = 1i/2;
        H1(k+1,k) = -1i/2;
        H3(k,k+1) = -1/2;
        H3(k+1,k) = -1/2;
        H4(k,k+1) = -1;
        H4(k+1,k) = -1;
    end
    k = i + (Ly-1)*Ly;
    k0 = i;
    H1(k,k0) = -1i/2;
    H1(k0,k) = 1i/2;
    H3(k,k0) = 1/2;
    H3(k0,k) = 1/2;
    H4(k,k0) = 1;
    H4(k0,k) = 1;
end

H = 3*(kron(H1,sigma_x) - kron(H2,sigma_y) + kron(H3,sigma_z)) + kron(H4,eye(2));

[phi,e] = eig(H);

phi_l1 = zeros(Ly,1);
phi_l2 = zeros(Ly,1);
phi_r1 = zeros(Ly,1);
phi_r2 = zeros(Ly,1);

for t = 1:Lx*Ly
    phi_l1 = phi_l1 + phi(1:Ly,t);
    phi_l2 = phi_l2 + phi(((Lx-1)*Ly+1):Lx*Ly,t);
    phi_r1 = phi_r1 + phi((Lx*Ly+1):((Lx+1)*Ly),t);
    phi_r2 = phi_r2 + phi(((2*Lx-1)*Ly+1):2*Lx*Ly,t);
end

phi_l1c = conj(phi_l1);
phi_l2c = conj(phi_l2);
phi_r1c = conj(phi_r1);
phi_r2c = conj(phi_r2);

J_l1 = zeros(1,Ly);
J_r1 = zeros(1,Ly);
J_l2 = zeros(1,Ly);
J_r2 = zeros(1,Ly);

for i = 1:Ly-1
    J_l1(i) = -1i*(phi_l1c(i)*phi_l1(i+1) - phi_l1(i)*phi_l1c(i+1));
    J_r1(i) = -1i*(phi_r1c(i)*phi_r1(i+1) - phi_r1(i)*phi_r1c(i+1));
    J_l2(i) = -1i*(phi_l2c(i)*phi_l2(i+1) - phi_l2(i)*phi_l2c(i+1));
    J_r2(i) = -1i*(phi_r2c(i)*phi_r2(i+1) - phi_r2(i)*phi_r2c(i+1));
end
y = 1:Ly-1;
plot(y,J_l1(1:Ly-1));
hold on
plot(y,J_r1(1:Ly-1));

toc;
