function [J,DJ,xstore] = calcJDJ(z,xb,h,ht,L,a,dt,dx,nt,nx,B,R,Q,ystore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Calculate full cost function and gradient
%  for 4D-Var with model error forcing, including background term
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

% Distribute contents of z over x_0 and eta_i's
xo = z(1:nx);

eta = zeros(nx,nt-1);
n=1; i=1;
while n<=(nt-1)
for k=(nx+1):(nx*nt) 
 eta(i,n) = z(k);
 i = i+1;
    if rem(k,nx) == 0
    n=n+1;
    i=1;
    end
end
end

xstore = zeros(nx,nt);  % State trajectory storage of assimilating model
xstore(:,1) = xo;

%-------------------------------------------------------------------------%
%% Calc cost function
J=0;

% Calc first term of cost func involving Jb
J=J+0.5*((xstore(:,1)-xb)'*B*(xstore(:,1)-xb)) ...
   +0.5*((ystore(:,1) - h*ht(1)*xstore(:,1))'*R*(ystore(:,1)-h*ht(1)*xstore(:,1)));

for n=2:nt

[mx]=adv1d(xstore(:,n-1),eta(:,n-1),L,a,1,dt,dx);

xstore(:,n) = mx;

J =J+0.5*((ystore(:,n)-h*ht(n)*xstore(:,n))'*R*(ystore(:,n)-h*ht(n)*xstore(:,n))) ...
    +0.5*(eta(:,n-1)'*Q*eta(:,n-1));
end

%-------------------------------------------------------------------------%

% Calc Gradient of cost function
djn = zeros(nx,nt);        % Gradient place holder
DJ = zeros(nx*nt,1);       % Gradient
lambda=zeros(nx,nt);       % Adjoint variable: each column is 
                           % 1 set of grid points at a given time

% First adjoint variable calculated outside loop
lambda(:,nt) = h*ht(nt)*R*(ystore(:,nt) - h*ht(nt)*xstore(:,nt));
djn(:,nt) = Q*eta(:,nt-1) - lambda(:,nt);

  for n=nt-1:-1:2      
   [mlambda]=adv1d_adj(lambda(:,n+1),eta(:,n),L,a,1,dt,dx);  
   % lambda used as initial cond's for adjoint model, to calc grad

   lambda(:,n) = mlambda + h*ht(n)*R*(ystore(:,n) - h*ht(n)*xstore(:,n));   

   djn(:,n) = Q*eta(:,n-1) - lambda(:,n);  
  end

% Final adjoint model run done seperately because \nabla J's
% first term of gradient involves background term, not model error term as above
  [mlambda]=adv1d_adj(lambda(:,2),eta(:,1),L,a,1,dt,dx);

  lambda(:,1) = mlambda + h*ht(1)*R*(ystore(:,1) - h*ht(1)*xstore(:,1));

  djn(:,1) = B*(xstore(:,1)-xb) - lambda(:,1);

% Allocate all values from gradient placeholder (djn) to DJ.
k=1;
    for n=1:nt
      for i=1:nx
          DJ(k) = djn(i,n);
          k=k+1;
      end
    end
end