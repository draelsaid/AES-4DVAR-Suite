function [M,D] = jord(a,d,small)
%Jord    Compute the Jordan decomposition of square matrix A
%        such that A*M = M*D.  The elements of A may be real
%        or complex floating point numbers.
%
%usage: [M, D] = jord(A, e, small)
%
%        The eigenvalues of A will appear on the diagonal of D
%        and the generalized eigenvectors of A are the columns of M.
%
%        The accuracy of this function is limited by the accuracy
%        in determining the exact eigenvalues of A.
%        For this reason this routine requires the user provide
%        the eigenvalues of A in the input vector, e.
%
%        This routine is used by the JORDAN function.
%        and uses the built in NULLD function.
%
%        Determining the Jordan form is a numerically ill conditioned
%        problem. If the condition number of A is large then you may
%        not get accurate results. Also,  if you have problems then
%        consider running the program again. A random combination of
%        Null Space vectors is used each time in an attempt to avoid
%        matrix singularity problems due to linearly dependent Null
%        Space subspace issues.
%        Also, consider changing the optional SMALL parameter.
%        Typically small ~= 2*sqrt(eps). If the program detects an error
%        then the program variables are saved in the JORD.MAT file for
%        inspection by the JordLook program attached below.
%
%
%tested under 5.3.1
%
%see also: sym/jordan eig hess schur

%Paul Godfrey
%pgodfrey@intersil.com
%Feb., 7, 2001

norma=sqrt(mean(mean(abs(a))));
tiny=norma*eps;

if nargin<3
   small=(1+norma)*sqrt(eps);
end

[r,c]=size(a);
if r~=c
   error('A matrix must be square!')
end
n=r;
I=eye(n);

if r==1
   D=a;
   M=1;
   return
end

if r<1
   D=[];
   M=[];
   return
end

condofa=cond(a);
if condofa>1e6
   Condition_number_of_A=condofa
   warning('Condition number of A is very large!');
end


d=d(:);
if size(d,1)~=n
   d=d
   error('Wrong number of eigenvalues provided!')
end

da=det(a);
dp=prod(d);
e=abs(abs(da)-abs(dp));
if e>sqrt(eps)
   disp(' ')
   warning('Eigenvalues are probably inaccurate!')
end

ds=flipud(sort(d));
sds=size(ds,1);

du=flipud(unique(ds));
sdu=size(du,1);

if sdu==sds
%  no repeated evalues at all
   [M,D]=eig(a);
   return
end

%  yes some repeated evalues
M=[];

for kk=1:sdu
       e=du(kk);
       ameig=sum(ismember(ds,e));% evalue alg mult
       a1=a-e*I;
    if ameig==1
       [u,s,v]=svd(a1);
       M=[M v(:,end)];
    else
       pp=0;
       ns=[];

           pp=pp+1;
           aa=I;

           for k=1:ameig
               aa=a1*aa;
               nn=size(nulld(aa,small),2);
               ns(k)=nn;
           end

       nsaa=[0; ns(:)]';
       dns=diff(nsaa);

       if max(dns)~=dns(1)
          Cond_of_A=cond(a)
          save jord
          M=I; D=I;
          error('Null space size error 1')
       end

%      determine the number and lengths of the eigen-chains
       clear ec;
       ec(1:dns(1))=1; % number of chains
       for k=2:length(dns)
           ec(1:dns(k))=ec(1:dns(k))+1;
       end

       if sum(ec)~=ameig
          Cond_of_A=cond(a)
          save jord
          M=I;D=I;
          error('Null space size error 2')
       end

%     ec(1) contains the length of the longest chain(s)
      k=1;
      
      clear jv;
      while k<= dns(1)
          p=find(ec==ec(k));
%         could be more than 1 chain of a given length
          if isempty(p)
             Cond_of_A=cond(a)
             save jord
             M=I; D=I;
             error('Null space size error 3');
          end
          aa=I;
          for m=1:ec(k)
             aa=aa*a1;
          end

%         [u,s,v]=svd(aa);
%         last cols of v are the nullspace

          pp=max(size(p));
          v=nulld(aa, small);         
          jv(:,p)=v*(rand(size(v,2),pp)-0.5)*16;
%         take a random(+-8) LC of the NSVs to make sure
%         that we use no repeats from v.--> Large cond(M)
%         This can happen since N1cN2cN3...  
          k=k+pp;
      end
      clear v;

      for k=1:dns(1) % do each chain
          v(:,1)=jv(:,k);% termination vector for kth chain
          for p=2:ec(k)
              v(:,p)=a1*v(:,p-1);
          end
          vv=fliplr(v(:,1:ec(k)));
          M=[M vv];
      end
   end
end

k=abs(det(M))^(-1/n);
M=k*M; % make det(m)=+-1
% if we did our job right then
% m is square and always non-singular
Mi=inv(M);
D=Mi*a*M;
d0=diag(D);
d1=diag(D,1);
D=diag(d0)+diag(d1,1);

return