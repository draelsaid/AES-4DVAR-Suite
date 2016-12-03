function [phi] = testgrad(cf,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Purpose: To test gradient of cost function
%
%  (c) 2012 University of Reading
%
%  Written by Adam El-Said
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(x,1); s = size(x,2);
z=randn(N,s);
dz = randn(N,s);

for i=1:16
alpha(i)=10^(1-i);      % Step sizes

[Jz,DJz] = cf(z);       % Store J(z) and J'(z)

zdz = z + alpha(i)*dz; 

[Jzdz,~] = cf(zdz);     % J(z+dz)

phi(i) = (Jzdz - Jz) / (alpha(i)*(reshape(dz,1,s*N)*DJz));
res(i) = abs( phi(i) - 1 );
end

% Plot variation of phi with the stepsize alpha
figure(11), loglog(alpha,abs(phi)); title('Verification of gradient calculation'), xlabel('\alpha'), ylabel('\Phi(\alpha)')

% Plot variation of residual with the stepsize alpha
figure(12), loglog(alpha,res), title('Variation of residual'), xlabel('\alpha'), ylabel('\Phi(\alpha)-1')
end