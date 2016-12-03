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
N = size(x,1); n = size(x,2);

if size(x,1)==N
  x = x(obsxo:obs_f_x:N,:);           % Spatial obs selection
else
end

if size(x,2)==n
  x = x(:,obsto:obs_f_t:n);           % Temporal obs selection
else
end

  Hx = x;                             % 4D H has now been applied to x 

end