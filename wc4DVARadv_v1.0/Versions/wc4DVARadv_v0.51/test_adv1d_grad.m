function [phi] = test_adv1d_grad(xb,h,ht,a,dt,dx,nt,nx,B,R,Q,ystore)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: To test gradient of cost function calcJDJ.m
%
%  (c) 2012  University of Reading
%
%  Written  by Adam El-Said
%
%  (*) This program requires alot of parameters, so its best to run when
%  w4dvar.m has been initialised properly. The call to this function will
%  be stored in the correct place (and commented out) in w4dvar.m

%  We have a perturbation such that:
%
%  dz = alpha * [ J' / ||J'|| ]  (*) This normalises the perturbation relative
%                                    to the function. Also 
%                                (*) alpha is the step size, [J' / ||J'||]
%                                    is the direction
%
%  [ J(z+dz) - J(z) ] / dz J'(z) = 1 + phi(alpha)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ic = 16;               % Counter

phi = zeros(ic,1);     % Function representing higher order terms of taylor expansion of J
alpha = zeros(ic,1);   % Step size
res = zeros(ic,1);     % [ 1 - phi(alpha) ]

z=0.1*randn((nx+1)*(nt+1),1);

for i=1:ic

% Calculate cost function and gradient with initial z
[J,DJ]=calcJDJ(z,xb,h,ht,a,dt,dx,nt,nx,B,R,Q,ystore);

alpha(i)=10^(1-i);      % Step sizes

Jz = J;                 % Store J(z)
DJz = DJ;               % Store J'(z)

nDJz = DJz/norm(DJz);   % Normalised gradient direction
deltaz = alpha(i)*nDJz; % Perturbation
z = z+deltaz;           % New z = z + dz

% Calculate cost function with new perturbed z
[J]=calcJDJ(z,xb,h,ht,a,dt,dx,nt,nx,B,R,Q,ystore);

Jzdeltaz = J;           % J(z+dz)

phi(i) = (Jzdeltaz - Jz) / (alpha(i)*nDJz'*DJz);
res(i) = abs( phi(i) - 1 );
end

%Plot variation of phi with the stepsizes alpha
figure(11)
semilogx(alpha,phi);
title('Verification of gradient calculation')
xlabel('\alpha')
ylabel('\Phi(\alpha)')

%Plot variation of residual with the stepsizes alpha
figure(12)
loglog(alpha,res)
title('Variation of residual')
xlabel('\alpha')
ylabel('\Phi(\alpha)-1')

end