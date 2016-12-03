function [M,D] = jordan(a)

%Jordan  Compute the Jordan decomposition of square matrix A
%        such that A*M = M*D.  The elements of A may be real
%        or complex floating point numbers.
%
%        The eigenvalues of A are on the diagonal of D
%        and the generalized eigenvectors of A are the columns of M.
%
%        This routine uses the Jord function.
%
%usage: [M, D] = jordan(A)
%
%tested under 5.3.1
%
%see also: sym/jordan mhelp('jordan') eig hess schur

%Paul Godfrey
%pgodfrey@intersil.com
%2-05-01

small=2*sqrt(eps);

[r,c]=size(a);
if r~=c
   error('A matrix must be square!')
end
n=r;

if n==1
   D=a;
   M=1;
   return
end

if n<1
   D=[];
   M=[];
   return
end

[m,d]=eig(hess(a));
d=sort(diag(d));
tiny=norm(a)*eps;

%zero out zero evalues
p=find(abs(d)<=tiny);
if ~isempty(p)
   d(p)=0;
end

%A*M=M*D

[M,D]=jord(a,d,small);

return