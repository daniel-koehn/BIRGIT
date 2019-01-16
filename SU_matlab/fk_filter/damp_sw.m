function [d] = damp_sw(d,offset,t,maxoffset_sw,vmin_sw,vmax_sw,gamma_sw)
% DAMP_SW: Damp surface waves in time domain
%
% Daniel Koehn
% Kiel, 05.09.2016

    [nt,nx]=size(d);

    for i=1:nx
         for j=1:nt

             % calculate velocity 
             v = abs(offset(i)/t(j));         

             if((v>=vmin_sw)&&(v<=vmax_sw)&&(abs(offset(i))<maxoffset_sw))
                d(j,i) = d(j,i) * exp(-gamma_sw);
             end

         end
    end
 
return;
