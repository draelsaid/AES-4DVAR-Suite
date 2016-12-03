function [phi] = testgrad(cf,x,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Purpose: To test gradient of cost function
%
%  (c) 2012 University of Reading
%
%  Written by Adam El-Said
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 2) 
N = size(x,1); s = size(x,2);
z=randn(N,s);
dz = randn(N,s);
else
N = size(x,1);
z=randn(N,s);
dz = randn(N,s);
end

for i=1:16
alpha(i)=10^(1-i);      % Step sizes

[Jz,DJz] = cf(z);       % Store J(z) and J'(z)

zdz = z + alpha(i)*dz; 

[Jzdz,~] = cf(zdz);     % J(z+dz)

phi(i) = (Jzdz - Jz) / (alpha(i)*(reshape(dz,1,s*N)*DJz));
res(i) = abs( phi(i) - 1 );
end

% Plot variation of phi with the stepsize alpha
figure, loglog(alpha,abs(phi)); title('Verification of gradient calculation'), xlabel('\alpha'), ylabel('\Phi(\alpha)')

% Plot variation of residual with the stepsize alpha
figure, loglog(alpha,res), title('Variation of residual'), xlabel('\alpha'), ylabel('\Phi(\alpha)-1')
end