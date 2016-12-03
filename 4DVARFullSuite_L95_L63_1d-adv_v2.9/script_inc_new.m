% Weak Constraint 4DVAR main sheet
% Adam El-Said, University of Reading (06/2013)

clear all
% Close any open figure windows
g=gcf;
close all

%% Graph position settings
ss=get(0,'ScreenSize'); fig_width=0.55*ss(3); fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height]; pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height]; pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3); pb=0.05*ss(4); pw=0.3*ss(3); ph=0.35*ss(4);
pos_4=[pl,pb,pw,ph]; pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph]; pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

N=40; n=10; dt=0.025; % 0.025 dt is a 3 hour time step (L95)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRUTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnlt = @l95; scht = @rk4nl;
xto = Mnl(rand(N,1),fnlt,scht,dt,n); % Model 'spin up' to obtain true initial state

[~,xttraj] = Mnl(xto,fnlt,scht,dt,n);
% xttraj now contains the TRUE trajectory.

% Model function handle and integration scheme setup
fnl = @l95;      schnl = @rk4nl;
ftl = @l95tl;    schtl = @rk4tl;
fadj= @l95adj;  schadj = @rk4adj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are now the propogators that will be used throughout the assimilation
fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 
fadjm = @(dy,y) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

% Covariance Matrices (TRUTH)
sigmabt = 0.1; Lsbt = 0.025;    % sigma and lengthscale for B
sigmaqt = 0.05; Lsqt = 0.01;    % sigma and lengthscale for Q

fcBt = @(N) soar(N,Lsbt,dt);            % B MATRIX function handle
fcQt = @(N) laplacian(N,Lsqt,dt);       % Q MATRIX function handle

% Background creation (model error backgrounds (\eta^b) are assumed 0 for now)
fDt = @(x) dmtx(x,sigmabt,sigmaqt,fcBt,fcQt);    % D MATRIX function handle
[~,pb] = fDt(xttraj); 
pb = pb + xttraj; 
pb(:,2:end)=0;   % Background x^b_0 is truth perturbed with noise

[~,xbtraj] = Mnl(pb(:,1),fnlt,scht,dt,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Observation creation parameters
sigmaot = 0.01;                  % True observation error variance

fRt = @(x) rmtx(x,sigmaot);      % Setting handle for use inside cost function
[~,y] = rmtx(xttraj,sigmaot);    % Observations are from the truth with some noise
y = y + xttraj;

obs_f_t = 0;      % Obs frequency (time) (Set to 0 if you want to observe ALL times)
obs_f_x = 1;      % Obs frequency (space) (Regular spacing)
obsxo = 1;        % obsxo determines observation starting point (space)
obsto = 1;        % obsto determines observation starting point (time)

fx = floor(N/obs_f_x); % Variable sizes determined by input

fH = @(y) Hop(y,obs_f_x,obs_f_t,obsxo,obsto);         % H operator
fHT = @(y) HopT(y,N,n+1,obs_f_x,obs_f_t,obsxo,obsto); % H operator Transposed

y = fH(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Covariance Matrices (ASSIMILATION)
sigmab = 0.1; Lsb = 0.025;    % sigmab and lengthscale for B
sigmaq = 0.05; Lsq = 0.01;    % sigmaq and lengthscale for Q
sigmao = 0.01;                % sigmao

fcB = @(N) soar(N,Lsb,dt);         % B MATRIX function handle
fcQ = @(N) laplacian(N,Lsq,dt);    % Q MATRIX function handle

fD = @(x) dmtx(x,sigmab,sigmaq,fcB,fcQ);
fR = @(x) rmtx(x,sigmao);

p = pb; x = pb;

% %% Gradient Test
%   cfdp = @(dp) wcdp(dp,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%   cfdx = @(dx) wcdx(dx,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%   phi = testgrad(cfdp,pb); 

%% Minimisation parameters
xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'100','1d-4'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr); 
max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

%% Polack-Ribiere Non-Linear Conjugate Gradient
s = size(pb,2);
p = reshape(p,s*N,1);
x = reshape(x,s*N,1);

[X_adp JXdp dJXdp itdp] = minimize_mod_crit(p,'wcdp',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);

X_adp = reshape(X_adp,N,s);
xatrajdp = Lopnlinv(X_adp,fnlm);

[X_adx JXdx dJXdx itdx] = minimize_mod_crit(x,'wcdx',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);

X_adx = reshape(X_adx,N,s);
xatrajdx = X_adx;

% Error calculations
% Relative norm error
for i=1:n+1
errordx(i) = norm(xttraj(:,i) - xatrajdx(:,i))/norm(xttraj(:,i));
end

for i=1:n+1
errordp(i) = norm(xttraj(:,i) - xatrajdp(:,i))/norm(xttraj(:,i));
end

% Relative error 
for i=1:N
    for j=1:n+1
      cerrordx(i,j) = abs(xttraj(i,j) - xatrajdx(i,j))/abs(xttraj(i,j));
    end
end

for i=1:N
    for j=1:n+1
      cerrordp(i,j) = abs(xttraj(i,j) - xatrajdp(i,j))/abs(xttraj(i,j));
    end
end

contourf(xatrajdx'), title('Contour plot Lorenz95 (dx Solution)'), ylabel('Time (n)'), xlabel('Variables (X_i)'); 
figure; contourf(xatrajdp'), title('Contour plot Lorenz95 (dp Solution)'), ylabel('Time (n)'), xlabel('Variables (X_i)');
figure; contourf(xttraj'), title('Contour plot Lorenz95 (TRUTH)'), ylabel('Time (n)'), xlabel('Variables (X_i)');
figure; contourf(xbtraj'), title('Contour plot Lorenz95 (Background)'), ylabel('Time (n)'), xlabel('Variables (X_i)'); 
figure; plot(errordp), title('Relative error (dp)'), xlabel('Time (n)'), ylabel('Error value')
figure; plot(errordx), title('Relative error (dx)'), xlabel('Time (n)'), ylabel('Error value')

figure; contourf(cerrordx), title('Contour plot Lorenz95 (dx errors)'), ylabel('Time (n)'), xlabel('Variables (X_i)'), colorbar;
figure; contourf(cerrordp), title('Contour plot Lorenz95 (dp errors)'), ylabel('Time (n)'), xlabel('Variables (X_i)'), colorbar;

%% Convergence plot dp
h3=figure('Position',pos_6); clf;
subplot(2,1,1)
semilogy(linspace(1,max(size(JXdp)),max(size(JXdp))),JXdp), title('Convergence of cost function (dp)'), xlabel('Iteration'), ylabel('Cost function')
subplot(2,1,2)
for i = 1 : itdp
  g =  dJXdp(:,i);
  graddp(i)=norm(g);
end
semilogy(linspace(1,max(size(graddp)),max(size(graddp))),graddp), title('Convergence of gradient (dp)'), xlabel('Iteration'), ylabel('Norm of gradient')

%% Convergence plot dx
h3=figure('Position',pos_6); clf; 
subplot(2,1,1)
semilogy(linspace(1,max(size(JXdx)),max(size(JXdx))),JXdx), title('Convergence of cost function (dx)'), xlabel('Iteration'), ylabel('Cost function')
subplot(2,1,2)
for i = 1 : itdx
  g =  dJXdx(:,i);
  graddx(i)=norm(g);
end
semilogy(linspace(1,max(size(graddx)),max(size(graddx))),graddx), title('Convergence of gradient (dx)'), xlabel('Iteration'), ylabel('Norm of gradient')
