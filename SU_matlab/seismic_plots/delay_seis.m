function [D_data] = delay_seis(delay,DT,ntr,nt,D_data,gamma_delay)
% DELAY_SEIS
% apply time delay to seismic data
% 
% Daniel Koehn
% Kiel, the 27th of September 2014

ndelay = round(delay./DT);
zdelay = zeros(ntr,ndelay);

% set values in delay zone to first sample of each trace
for i=1:ntr
   zdelay(i,:) = D_data(1,i);
end

D_data = [zdelay' ; D_data];

tmp = D_data(1:nt,:);
clear D_data;
D_data = tmp;
clear tmp;

% damp data before begin of original dataset 
for i=1:ntr
   for j=1:ndelay
       D_data(j,i)=D_data(j,i).*exp(-gamma_delay.*(j-ndelay).^2);
   end
end

clear zdelay;

end

