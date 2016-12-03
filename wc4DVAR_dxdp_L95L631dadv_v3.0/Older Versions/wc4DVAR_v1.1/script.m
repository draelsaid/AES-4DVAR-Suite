clear all
close all

N=40;n=20;dt=0.025;
x=randn(N,n);

fnlt = @l95;
scht = @rk4nl;
xto = Mnl(rand(N,1),fnlt,scht,dt,n); % Model 'spin up' to obtain true initial state

[~,xttraj] = Mnl(xto,fnlt,scht,dt,n);
% xttraj now contains the TRUE trajectory.

% Model function handle setup
fnl = @l95; schnl = @rk4nl;
ftl = @l95tl; schtl = @rk4tl;
fadj = @l95adj; schadj = @rk4adj;

fnlm = @(x) Mnl(x,fnl,schnl,dt,n);
ftlm = @(x,dx) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n); 
fadjm = @(y,dy) Madj(dy,y,fadj,fnl,schadj,schnl,dt,n);

% Covariance Matrices (assim)
sigmab = 1; Lsb = 0.1;
sigmaq = 1; Lsq = 0.2;

fcB = @(N) soar(N,Lsb,1);
fcQ = @(N) soar(N,Lsb,1);

% Background creation (model error backgrounds are assumed 0 for now)
fD = @(x) dmtx(x,sigmab,sigmaq,fcB,fcB); % Setting handle for use inside cost function
[~,pb] = fD(xttraj); pb(:,2:end)=0;      % Background x^b_0 is truth perturbed with some noise

[~,xbtraj] = fnlm(pb(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Observation creation parameters
sigmar = 1;

fR = @(x) rmtx(x,sigmar);       % Setting handle for use inside cost function
[~,y] = rmtx(xttraj,sigmar);    % Observations are from the truth with some noise

obs_f_t = 2;      % Obs frequency (time) (Set to 0 if you want to observe ALL times)
obs_f_x = 2;      % Obs frequency (space)
obsxo = 1;        % obsxo determines where you want to start observing from (space)
obsto = 1;        % obsto determines where you want to start observing from (time)

fx = floor(N/obs_f_x); % Variable sizes determined by input

fH  = @(y) Hop(y,obs_f_x,obs_f_t,obsxo,obsto);        % H operator
fHT = @(y) HopT(y,N,n+1,obs_f_x,obs_f_t,obsxo,obsto); % H operator transposed

y = fH(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=pb;

cf = @(z) wcdp(z,pb,y,n,fH,fHT,fD,fR,fnl,fadj,schnl,schadj,dt);

phi = testgrad(cf,pb); % Gradient Test

%% Minimisation parameters
xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'30','1d-6'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr); 
max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

%% Polack-Ribiere Non-Linear Conjugate Gradient

 [X_a JX dJX it] = minimize_mod_crit(p,cf,max_iterations,tolerance);
 
 
 