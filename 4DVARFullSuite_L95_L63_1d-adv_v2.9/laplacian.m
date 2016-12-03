function [Cl,gamma] = laplacian(N,Ls,dx)

q = (Ls^4) / (2*dx^4); % Constant dependent on lengthscale and gridspacing

% First row of Laplacian
r = zeros(1,N); r(1) = 1+6*q; r(2) = -4*q; r(3) = q; r(N-1) = q; r(N) = -4*q;
Cl = toeplitz(r); % Create Laplacian inverse using the first row

Clinv = inv(Cl);  % Constant to ensure maximum element of first row doesn't exceed 1
gamma = Clinv(1,1);
% Produces same matrix in a different way
%r2 = zeros(1,N); r2(1) = -2; r2(2) = 1; r2(N) = 1;
%Cl2 = gamma*(eye(N,N) + q*(toeplitz(r2))^2);
Cl = gamma^(-1)*Clinv;
end
