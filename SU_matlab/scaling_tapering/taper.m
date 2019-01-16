function [dout,win] = taper(d,perc_beg,perc_end);
%TAPER: Apply a triangular  taper to the beg/end of traces.
% 
%  IN   d(nt,nx):      data (columns are traces)
%       perc_beg:      percetage of data to be tapered at the beggining
%       perc_end:      "                                      end
%
%  OUT  dout(nt,nx):   data with beg/end tapered
%       win(nt):       the taper
%
%  Example
%
%    d = make_traces; wigb([d,taper(d,50,50)]);
%
%  Copyright (C) 2008,  Signal Analysis and Imaging Group
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


 [nt,nx] = size(d);

 i1 = floor(nt*perc_beg/100)+1;
 i2 = floor(nt*perc_end/100)+1;
   
 win = [ [1:1:i1]/i1,ones(1,nt-i1-i2),[i2:-1:1]/i2];

 for k=1:nx
  dout(:,k) = d(:,k).*win';
 end;


