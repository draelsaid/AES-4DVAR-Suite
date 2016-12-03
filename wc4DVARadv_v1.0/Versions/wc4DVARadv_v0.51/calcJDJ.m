function [DJ] = calcJDJ(z,xb,h,ht,a,dt,dx,nt,nx,B,R,Q,ystore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate full 4DVAR cost function and gradient
%  for with model error forcing
%
%  Written by Adam El-Said
%
%  List of main variables
%
%% Assimilation variables
%    xb:         Background term
%    h:          Spatial component of observation operator
%    ht:         Assertion of observation times
%    B:          Background weighting matrix (inversed)
%    R:          Observation weighting matrix (inversed)
%    Q:          Model Error weighting matrix (inversed)
%    [z]:        Current iterate
%    [ystore]:   Observation values
%
%% Model variables
%    L:          Physical domain size
%    a:          Velocity in advection equation
%    dt:         Time step size
%    dx:         Spatial step size
%    nt:         Number of time steps
%    nx:         Number of gridpoints
%
%% Output:
%    [J,DJ]: Cost function and gradient
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% J_b, J_o and J_q Term Control %% [ Jq = 0 => STRONG CONSTRAINT 4DVAR ]
Jb = 1;
Jo = 1;
Jq = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distribute contents of z over x_0 and eta_i's
Z = reshape(z,nx+1,nt+1);
xo = Z(:,1);
eta = Z(:,2:nt+1);

xstore = zeros(nx+1,nt+1);  % State trajectory storage of assimilating model
xstore(:,1) = xo;

%-------------------------------------------------------------------------%
%% Calc cost function
J=0;

% Calc first term of cost func involving Jb
J=J+Jb*[0.5*((xstore(:,1)-xb)'*B*(xstore(:,1)-xb))] ...
   +Jo*[0.5*((ystore(:,1) - h*ht(1)*xstore(:,1))'*R*(ystore(:,1)-h*ht(1)*xstore(:,1)))];

for n=2:nt+1

[mx]=adv1d(xstore(:,n-1),eta(:,n-1),nx,a,1,dt,dx);

xstore(:,n) = mx;

J =J+Jo*[0.5*((ystore(:,n)-h*ht(n)*xstore(:,n))'*R*(ystore(:,n)-h*ht(n)*xstore(:,n)))] ...
    +Jq*[0.5*(eta(:,n-1)'*Q*eta(:,n-1))];
end

%-------------------------------------------------------------------------%

% Calc Gradient of cost function
djn = zeros(nx+1,nt+1);        % Gradient place holder
DJ = zeros((nx+1)*(nt+1),1);   % Gradient
lambda=zeros((nx+1),(nt+1));   % Adjoint variable

% First adjoint variable calculated outside loop
lambda(:,nt+1) = Jo*[h*ht(nt+1)*R*(ystore(:,nt+1) - h*ht(nt+1)*xstore(:,nt+1))];
djn(:,nt+1) = Jq*[(Q*eta(:,nt) - lambda(:,nt+1))];

  for n=nt:-1:2
   [mlambda]=adv1d_adj(lambda(:,n+1),eta(:,n),nx,a,1,dt,dx);
   % lambda used as initial cond's for adjoint model, to calc grad

   lambda(:,n) = mlambda + Jo*[h*ht(n)*R*(ystore(:,n) - h*ht(n)*xstore(:,n))];   

   djn(:,n) = Jq*[(Q*eta(:,n-1) - lambda(:,n))];  
  end

% Final adjoint model run done seperately because \nabla J's
% first term of gradient involves background term, not model error term as above
  [mlambda]=adv1d_adj(lambda(:,2),eta(:,1),nx,a,1,dt,dx);

  lambda(:,1) = mlambda + Jo*[h*ht(1)*R*(ystore(:,1) - h*ht(1)*xstore(:,1))];

  djn(:,1) = Jb*[B*(xstore(:,1)-xb) - lambda(:,1)];
  
% Allocate all values from gradient placeholder (djn) to DJ.
DJ = reshape(djn,(nt+1)*(nx+1),1);

% % For solvers that require the form Ax = b
% b = zeros((nx+1)*(nt+1),1);
% b(1:nx+1) = Jb*B*xb+lambda(:,1);
% b(nx+2:(nx+1)*(nt+1)) = reshape(lambda(:,2:nt+1),nt*(nx+1),1);
% 
% Ax = DJ - b;
end