function [dout] = gain(d,dt,option1,parameters,option2);
%GAIN: Gain a group of traces.
%
%  [dout] = gain(d,dt,option1,parameters,option2);
%
%  IN   d(nt,nx):   traces
%       dt:         sampling interval
%       option1 = 'time' parameters = [a,b],  gain = t.^a . * exp(-bt)
%               = 'agc' parameters = [agc_gate], length of the agc gate in secs
%       option2 = 0  No normalization
%               = 1  Normalize each trace by amplitude
%               = 2  Normalize each trace by rms value
%
%  OUT  dout(nt,nx): traces after application of gain function
%
%
%  Example:
%
%    d = hyperbolic_events; dout = gain(d,0.004,'agc',0.05,1);
%    wigb([d,dout]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%     


 [nt,nx] = size(d);

 if strcmp(option1,'time');   % Geometrical spreading-like gain

  a = parameters(1);
  b = parameters(2);
  t = (0:1:nt-1)*dt; 
  tgain = (t.^a).*exp(b.*t);

  for k = 1:nx
   dout(:,k)  = d(:,k).*tgain';
  end;

 end

 if strcmp(option1,'agc');    % AGC 

  L = parameters(1)/dt+1;
  L = floor(L/2);
  h = triang(2*L+1);

  for k = 1:nx
   aux =  d(:,k);
   e = aux.^2;
   rms = sqrt(conv2(e,h,'same'));
   epsi = 1.e-10*max(rms);
   op = rms./(rms.^2+epsi);
   dout(:,k) = d(:,k).*op;
   end
  end

 if option2==1;                % Normalize by amplitude 

   for k = 1:nx
    aux =  dout(:,k);
    amax = max(abs(aux));
    dout(:,k) = dout(:,k)/amax;
   end

 end


 if option2==2;                % Normalize by rms 

   for k = 1:nx
    aux =  dout(:,k);
    amax =  sqrt(sum(aux.^2)/nt);
    dout(:,k) = dout(:,k)/amax;
   end

 end

