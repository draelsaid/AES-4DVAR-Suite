% Weak Constraint 4DVAR main sheet
% Adam El-Said, University of Reading (06/2013)

clear all
% Close any open figure windows
g=gcf;
close all

%% Graph position settings
ss=get(0,'ScreenSize'); fig_width=0.55*ss(3); fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height]; pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height]; pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3); pbb=0.05*ss(4); pw=0.3*ss(3); ph=0.35*ss(4);
pos_4=[pl,pbb,pw,ph]; pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph]; pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

% Model settings
N=100; n=100; % Size of spatial and temporal domains
dt = 0.005;   % 0.025 dt => 3 hour time step (L95), 1 time unit is 5 days
dx = 1/N;     % Spatial interval size

%% TRUTH

fnlt = @adv1d; scht = @pde;     % Non-lin model and scheme func handle
% True Covariance Matrices D and R
sigmabt = 0.25; Lsbt = 0.05;    % True sigma_b and lengthscale for C_B
sigmaqt = 0.125; Lsqt = 0.025;  % True sigma_q and lengthscale for C_Q
sigmaot = 0.125;                % True sigma_o
fcBt = @(N) soar(N,Lsbt,dx);                  % B MATRIX function handle
fcQt = @(N) laplacian(N,Lsqt,dx);             % Q MATRIX function handle
fDt = @(x) dmtx(x,sigmabt,sigmaqt,fcBt,fcQt); % D MATRIX func handle
fRt = @(x) rmtx(x,sigmaot);                   % R MATRIX function handle

%% Initial Conditions

% Gaussian-shaped curve (for 1D advection equation ONLY)
xc = (N*dx)/2; sigma = 0.2; Ampl = 2; % Mode, Standard dev and Amplitude
xs = linspace(0,1,N);                 % Physical space index 

for i=1:N
    xto(i) = Ampl*exp(-(xs(i)-xc).^2/sigma^2);
end

% % For chaotic ODE models
% xto = Mnl(rand(N,1),fnlt,scht,dt,50); % Model 'spin up' to obtain true initial state

% Truth generation
[~,etat] =  fDt(rand(N,n+1));       % True model errors
bias = 0.1; etat = bias + etat;       % Model error bias
etat(:,1) = 0;                      % Model error time adj(eta_0 not exist)
[~,xt] = Mnl(xto,fnlt,scht,dt,n);   % Truth model trajectory

xttraj = xt + etat;

% % For advection equation; solution run through time
% for i=1:n
%  figure(10), plot(xs,xttraj(:,i),'b'), title('Truth run'), xlabel('Position'), legend('Analysis from t_0 to t_{now}'), ylabel('Value')
%  drawnow
% end

%% ASSIMILATION
% Model function handle and integration scheme setup
fnl = @adv1d;      schnl = @pde;
ftl = @adv1d;      schtl = @pde;
fadj= @adv1d_adj; schadj = @pde;

% Assimilation model propagators 
fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 
fadjm = @(dy,y) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

% Background creation [truth plus some noise ~N(0,D)]
[~,xttrrajerror] = fDt(xttraj);  
pb = xttrrajerror + xttraj;      
pb(:,2:end)=0;                   % Background model error eta^b=0

[~,xbtraj] = Mnl(pb(:,1),fnlt,scht,dt,n);

%% Observation creation parameters

[~,y] = fRt(xttraj);    % Obs = truth + noise ~N(0,R)
y = y + xttraj;

obs_f_t = 5;      % Obs frequency (time) (Set=0 to observe ALL times)
obs_f_x = 5;      % Obs frequency (space) (Regular spacing)
obsxo = 1;        % Determines observation starting point (space)
obsto = 1;        % Determines observation starting point (time)

fx = floor(N/obs_f_x); % Variable sizes determined by input

fH = @(y) Hop(y,obs_f_x,obs_f_t,obsxo,obsto);   % H operator
fHT = @(y) HopT(y,obs_f_x,obs_f_t,obsxo,obsto); % H Transposed

y = fH(y); % H operator selects Obs, both spatially and temporally

%% Covariance Matrices (ASSIMILATION)

sigmab = 0.25; Lsb = 0.05;     % Estimated sigmab and lengthscale for B
sigmaq = 0.125; Lsq = 0.025;   % Estimated sigmaq and lengthscale for Q
sigmao = 0.125;               % Estimated sigmao

fcB = @(N) soar(N,Lsb,dx);         % C_B MATRIX function handle
fcQ = @(N) laplacian(N,Lsq,dx);    % C_Q MATRIX function handle

fD = @(x) dmtx(x,sigmab,sigmaq,fcB,fcQ); % D MATRIX function handle
fR = @(x) rmtx(x,sigmao);                % R MATRIX function handle

p = pb; x = pb;

% %% Gradient Test
%   cfdp = @(dp) wcdp(dp,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%   cfdx = @(dx) wcdx(dx,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%   phi = testgrad(cfdp,pb); 

%% Minimisation

% Minimisation parameters
xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'100','1d-4'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr); 
max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

s = size(pb,2);
p = reshape(p,s*N,1);
x = reshape(x,s*N,1);

% Minimisation procedure for dp
[X_adp JXdp dJXdp itdp] = minimize_mod_crit(p,'wcdp',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);

X_adp = reshape(X_adp,N,s);
xatrajdp = Lopnlinv(X_adp,fnlm);

% Minimisation procedure for dx
[X_adx JXdx dJXdx itdx] = minimize_mod_crit(x,'wcdx',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);

X_adx = reshape(X_adx,N,s);
xatrajdx = X_adx;

%% Trajectory time adjustment
xttraj = xttraj(1:N,2:end);
xbtraj = xbtraj(1:N,2:end);
xatrajdx = xatrajdx(1:N,2:end);
xatrajdp = xatrajdp(1:N,2:end);

%% Error calculations

% Relative total error
terrordx = norm(xttraj - xatrajdx) / norm(xttraj)
terrordp = norm(xttraj - xatrajdp) / norm(xttraj)

% Relative error of each time interval
for i=1:n
errordx(i) = norm(xttraj(:,i) - xatrajdx(:,i))/norm(xttraj(:,i));
end

for i=1:n
errordp(i) = norm(xttraj(:,i) - xatrajdp(:,i))/norm(xttraj(:,i));
end

% Relative absolute error spatially and temporally
for i=1:N
    for j=1:n
      cerrordx(i,j) = abs(xttraj(i,j) - xatrajdx(i,j))/abs(xttraj(i,j));
    end
end

for i=1:N
    for j=1:n
      cerrordp(i,j) = abs(xttraj(i,j) - xatrajdp(i,j))/abs(xttraj(i,j));
    end
end

%% L95 Contour plots
contourf(xatrajdx'), title('Contour plot Lorenz95 (J(x) Solution)'), ylabel('Time (n)'), xlabel('Variables (X_i)'); 
figure; contourf(xatrajdp'), title('Contour plot Lorenz95 (J(p) Solution)'), ylabel('Time (n)'), xlabel('Variables (X_i)');
figure; contourf(xttraj'), title('Contour plot Lorenz95 (TRUTH)'), ylabel('Time (n)'), xlabel('Variables (X_i)');
figure; contourf(xbtraj'), title('Contour plot Lorenz95 (Background)'), ylabel('Time (n)'), xlabel('Variables (X_i)'); 

%% J(x), J(p) spatio-temporal error plots
figure; 
e1 = subplot(2,1,1);
contourf(cerrordx), title('Contour plot Lorenz95 (J(x) errors)'), xlabel('Time (n)'), ylabel('Variables (X_i)'), colorbar;
e2 = subplot(2,1,2);
contourf(cerrordp), title('Contour plot Lorenz95 (J(p) errors)'), xlabel('Time (n)'), ylabel('Variables (X_i)'), colorbar;
%set([e1 e2],'clim',[0 max( max(max(log(cerrordp))),  max(max(log(cerrordx))) )]);

%% J(x), J(p) relative space-averaged error plots vs time 
figure;
plot(errordp), axis([-inf,inf,-inf,inf]), title('Relative error'), xlabel('Time (n)'), ylabel('Error value')
hold; plot(errordx,'r'), legend('Errors J(p) formulation', 'Errors J(x) formulation'); hold;

%% Convergence plot dp
h3=figure('Position',pos_6); clf;
subplot(2,1,1)
semilogy(linspace(1,max(size(JXdp)),max(size(JXdp))),JXdp), title('Convergence of cost function, J(p)'), xlabel('Iteration'), ylabel('Cost function'), axis([1,inf,0,inf]),
subplot(2,1,2)
for i = 1 : itdp
  g =  dJXdp(:,i);
  graddp(i)=norm(g);
end
semilogy(linspace(1,max(size(graddp)),max(size(graddp))),graddp), title('Convergence of gradient, J(p)'), xlabel('Iteration'), ylabel('Norm of gradient'), axis([1,inf,0,inf]),

%% Convergence plot dx
h3=figure('Position',pos_6); clf; 
subplot(2,1,1)
semilogy(linspace(1,max(size(JXdx)),max(size(JXdx))),JXdx,'r'), title('Convergence of cost function, J(x)'), xlabel('Iteration'), ylabel('Cost function'), axis([1,inf,0,inf]),
subplot(2,1,2)
for i = 1 : itdx
  g =  dJXdx(:,i);
  graddx(i)=norm(g);
end
semilogy(linspace(1,max(size(graddx)),max(size(graddx))),graddx,'r'), title('Convergence of gradient, J(x)'), xlabel('Iteration'), ylabel('Norm of gradient'), axis([1,inf,0,inf]),
