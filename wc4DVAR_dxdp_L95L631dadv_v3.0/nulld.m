function Z = nulld(A, small)
%NULLD   Null space.
%   Z = NULL(A) is an orthonormal basis for the null space of A obtained
%   from the singular value decomposition.  That is,  A*Z has negligible
%   elements, size(Z,2) is the nullity of A, and Z'*Z = I.
%
%   Uses a modified singular value tolerance formula and method.
%
%
%   See also SVD, ORTH, RANK, RREF.

 
%Paul Godfrey
%pgodfrey@intersil.com
%Feb. 7, 2001

% small is ~ 2*sqrt(eps);
% used to change zero values into small*nonzero_min()
% and used to detect all zero singular values <small


norma=sqrt(mean(mean(abs(A))));
tiny=norma*eps;

if nargin<2
   small=(1+norma)*sqrt(eps);
end

[m,n] = size(A);
if m~=n
   error('Matrix must be square!')
end

p=find(abs(A)<tiny);
if ~isempty(p)
   A(p)=0;
end

[U,S,V] = svd(A,0);
S=diag(S);
s=S;
norma=max(s);

% normalize s
smax=max(s);
if smax==0;smax=1;end
s=s/smax;
snorm=s;

% if s is all zeros
t=find(s>0);
if isempty(t);Z=V;return;end

p=find(s<tiny);
if ~isempty(p)
   s(p)=tiny;
end

%duplicate smallest non-zero value*small
%so log(s) won't goof up
%not needed because of above code
p=find(s==0);
if ~isempty(p)
   s(p)=small*min(s(t));
end

% find the largest change between consecutive s values
% this should be the spot where the zero singular values begin.
logs=log10(s);
sdifflog=-diff(logs);
%good place to break for debug
smax=max(sdifflog);
r=find(sdifflog==smax);

% use the last index
% or 0 for Inf or NaN
if min(size(r))>0
   r=r(end);
else
   r=0;
end

% the above method can still fail for matrices with
% uniformly decreasing s values but it's more consistent
% than the Matlab null command.

% the null space vectors are found
% in the last few cols of V
Z = V(:,r+1:n);

% Fix degenerate cases to match Matlab's Null output

% Scaled Identity matrix has no null space at all
if snorm==ones(n,1);  Z=zeros(n,0); end

% All zero matrix has V as it's null space
%Is this really an all zero vector? 
if max(S)<=small; Z=V; end

% save data for debugging in the nulld.mat file
%save nulld

return


%JordLook  JORD.MAT data analyses program
%
%see also: Jordan, Jord

%Paul Godfrey
%pgodfrey@intersil.com
%Feb. 7, 2001

clc
clear all
load jord
 
aa=I;
for k=1:n
    aa=aa*a1;
    nsv=nulld(aa);
    sk(1,k)=size(nsv,2);
    sk(2:(n+1),k)=svd(aa);
end

disp('Null Space size and singular values') 
Data=sk
 
return