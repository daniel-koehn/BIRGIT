function [s] = stackgather(d,N);
%STACKGATHER: a program to stack (with normalization) one gather.
%             s(t) = [sum_x d(x,t) ]/N(t);
%
%  [s] = stackgather(d,N);
%
%  IN   d(nt,nh): data (gather)
%       N(nt,1):  normalization factor for each time 
%
%  OUT  s:        stack (a normalized spatial average of traces) 
%
%  If N is not provided, we divide by number of traces
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



 [nt,nh] = size(d);

 M = nargin;

 if M<2; N = nh*ones(nt,1); end;

 s = sum(d,2);

 for k = 1:nt;

   if N(k,1)>0; s(k,1) = s(k,1)/N(k,1); else;
                s(k,1) = 0.; end;

 end;
