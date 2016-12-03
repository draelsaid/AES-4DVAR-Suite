f = figure('Position',[440 500 461 146]);

% create the data
%d(:,1) = storeitdp;
%d(:,2) = storeitPdp;
d(:,1) = storeitdx;

%d(:,3) = storeerrdp;
d(:,2) = storeerrdx;
%d(:,4) = storeerrPdp;

%d(:,5) = storekappaSp;
%d(:,6) = storekappaPSp;
d(:,3) = storekappaSx;
%d(:,7) = storekappaD;
%d(:,8) = storehesssize;

d(:,4) = sigmabsigmaq;
d(:,5) = sigmabsigmao;
d(:,6) = sigmaqsigmao;

d(all(d==0,2),:)=[]
% % Create the column and row names in cell arrays 
% cnames = {'Jpit','Jperr','Jxit','Jxerr','sbsq','sbso','sqso'};
% 
% % Create the uitable
% t = uitable(f,'Data',d,...
%             'ColumnName',cnames);