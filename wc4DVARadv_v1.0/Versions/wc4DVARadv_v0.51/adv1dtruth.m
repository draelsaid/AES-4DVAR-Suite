function [mx]=adv1dtruth(xo,nx,a,tstep,dt,dx,etat)

%% 1D Linear Advection Model (PERFECT MODEL)
%  solved by an upwind scheme with periodic boundary condition

Tnew = zeros(nx+1,1);  % Array for next time step
T=xo;                  % Set initial conditions
C = (a*dt)/(dx);       % Courant number

%% Time Loop
for n=1:tstep

% Periodic Boundary conditions
if a>0
    Tnew(1) = (1-C)*T(1) + C*T(nx+1) + etat(1)*dt;
else
    Tnew(nx+1) = (1+C)*T(nx+1) + (-C)*T(1) + etat(nx)*dt;
end

% Upwind Scheme
if a>0
    for i=2:nx+1
        Tnew(i) = C*T(i-1) + (1-C)*T(i) + etat(i)*dt;
    end
 else
     for i=1:nx
         Tnew(i) = (1+C)*T(i) + (-C)*T(i+1) + etat(i)*dt;
     end
 end

% Update
    for i=1:nx+1
     T(i) = Tnew(i);
    end

end

mx=T;

end
