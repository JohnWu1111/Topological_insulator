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

phi_xy1 = zeros(Lx,Ly);
phi_xy2 = zeros(Lx,Ly);
Jx1 = zeros(Lx,Ly);
Jx2 = zeros(Lx,Ly);
Jy1 = zeros(Lx,Ly);
Jy2 = zeros(Lx,Ly);
for t = 1:Lx*Ly
    for i = 1:Lx
        for j = 1:Ly
            phi_xy1(i,j) = phi((i-1)*Lx+Ly,t);
            phi_xy2(i,j) = phi(Lx*Ly+(i-1)*Lx+Ly,t);
        end
    end
    
    phi_xy1c = conj(phi_xy1);
    phi_xy2c = conj(phi_xy2);
    
    for i = 1:Lx-1
        for j = 1:Ly
            Jx1(i,j) = Jx1(i,j)-1i*(phi_xy1c(i,j)*phi_xy1(i+1,j) - phi_xy1(i,j)*phi_xy1c(i+1,j));
            Jx2(i,j) = Jx2(i,j)-1i*(phi_xy2c(i,j)*phi_xy2(i+1,j) - phi_xy2(i,j)*phi_xy2c(i+1,j));
        end
    end
    
    for i = 1:Lx
        for j = 1:Ly-1
            Jy1(i,j) = Jy1(i,j)-1i*(phi_xy1c(i,j+1)*phi_xy1(i,j) - phi_xy1(i,j)*phi_xy1c(i,j+1));
            Jy2(i,j) = Jy2(i,j)-1i*(phi_xy2c(i,j+1)*phi_xy2(i,j) - phi_xy2(i,j)*phi_xy2c(i,j+1));
        end
    end
end

toc;
