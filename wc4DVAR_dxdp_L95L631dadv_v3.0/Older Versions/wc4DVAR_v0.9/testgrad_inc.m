function [phi] = testgrad_inc(cf,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Purpose: To test gradient of cost function (incremental)
%
%  (c) 2012 University of Reading
%
%  Written by Adam El-Said
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(x,1); s = size(x,2);
dz=randn(N,s); ddz = randn(N,s); 

for i=1:16
alpha(i)=10^(1-i);     

[Jdz,DJdz] = cf(dz,ddz);      

dzddz = dz + alpha(i)*ddz; 

[Jdzddz,~] = cf(dzddz,ddz);    

phi(i) = (Jdzddz - Jdz) / (alpha(i)*(reshape(ddz,1,s*N)*DJdz));
res(i) = abs( phi(i) - 1 );
end

% Plot variation of phi with the stepsize alpha
figure(11), loglog(alpha,abs(phi)); title('Verification of gradient calculation'), xlabel('\alpha'), ylabel('\Phi(\alpha)')

% Plot variation of residual with the stepsize alpha
figure(12), loglog(alpha,res), title('Variation of residual'), xlabel('\alpha'), ylabel('\Phi(\alpha)-1')
end