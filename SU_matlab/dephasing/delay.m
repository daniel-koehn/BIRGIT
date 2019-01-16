function dout = delay(d1,d2,max_delay);
%DELAY: Delay d2 toward d1
%
%  [dout] = delay(d1,d2,max_delay);
%
%  IN   d1:        Reference data (vector or matrix, traces are columns)
%       d2:        Data to delay (d1 and d2 must have the same size)
%       max_delay: Max x-correlation lag
%
%  OUT  dout:      Delayed version of d2 toward reference data d1
%
%  Note: Use this program to apply a +/- time shift after deconvolution.
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

  [nt,nx] = size(d1);

  r = zeros(2*max_delay+1,1);

   for k = 1:nx
    temp = xcorr(d1(:,k),d2(:,k),max_delay);
    r = r + temp;
   end;

 [rmax,Lag] = max(r);
 time_to_move = Lag-max_delay-1

 if time_to_move>0;
  dout = [zeros(time_to_move,nx);d2];
  dout = dout(1:nt,:);
 end;
    
 if time_to_move ==0; dout=d2; end;

 if time_to_move<0; 
  ll = -time_to_move+1;
  dout = [d2(ll:nt,:);zeros(ll-1,nx)];
 end;
 
 return;
