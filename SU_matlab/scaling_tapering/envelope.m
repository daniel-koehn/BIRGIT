function [e] = envelope(d);
%ENVELOPE: The envelope of a group of traces.
%
%  e = envelope(d);
%
%  IN   d(nt,np):  traces
%
%  OUT  e(nt,np):  envelopes
%
%
%  Example
%
%    [d,r,t] = make_traces; d=d(:,1); e = envelope(d); 
%    plot(t,d,t,e);
%
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



[nt,np] = size(d);

  for k=1:np

  aux = d(:,k);
  u = hilbert(aux);
  e(:,k) = abs(u);

 end;



