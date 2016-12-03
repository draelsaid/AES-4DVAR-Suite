function Hx = HopT(x,N,n,obs_f_x,obs_f_t,obsxo,obsto)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H Operator Transposed                                           %
%     [Hx] - Output variables observed in space and time          %
%        x - Matrix containing observations (Col-Time, Row-Space) % 
%  obs_f_x - Frequency of observations in space                   %
%  obs_f_t - Frequency of observations in time                    %
%    obsxo - Starting point of spatial obs                        %   
%    obsto - Starting point of temporal obs                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx = floor(N/obs_f_x);

H = zeros(fx,N);

j=obsxo;
for  i = 1:fx
   H(i,j) = 1;
   j = j + obs_f_x;
end

  x = H'*x;
  x(:,obsto:obs_f_t:n) = [];
  Hx = x;

end