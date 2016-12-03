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
obst = size(obsto:obs_f_t:n); ft = obst(2); % Determines number of temporal obs
obsx = size(obsxo:obs_f_x:N); fx = obsx(2); % Determines number of spatial obs
H = zeros(fx,N); h = zeros(ft,n); 

j=obsxo;
for  i = 1:fx
   H(i,j) = 1;
   j = j + obs_f_x;
end

x=H'*x;                % Spatial H

j=obsto;
for  i = 1:ft
   h(i,j) = 1;
   j = j + obs_f_t;
end

 x = x*h;              % Temporal h
Hx = x;                % 4D H^T has now been applied to x 

end