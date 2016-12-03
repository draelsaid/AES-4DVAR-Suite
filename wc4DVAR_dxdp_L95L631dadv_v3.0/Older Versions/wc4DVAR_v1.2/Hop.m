function Hx = Hop(x,obs_f_x,obs_f_t,obsxo,obsto)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H Operator                                                      %
%     [Hx] - Output variables observed in space and time          %
%        x - Matrix containing observations (Col-Time, Row-Space) % 
%  obs_f_x - Frequency of observations in space                   %
%  obs_f_t - Frequency of observations in time                    %
%    obsxo - Starting point of spatial obs                        %   
%    obsto - Starting point of temporal obs                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(x,1); n = size(x,2); % Variable sizes (space and time) determined by input matrix

j = 1;
for i = obsxo:obs_f_x:N
    H(i,obs_f_x*j) = 1.;
    j = j + 1;
end

  x = H*x; 
  x(:,obsto:obs_f_t:n+1) = 0;
  Hx = x;
end