function sgray(alpha)
%SGRAY: Non-linear transformation of a gray colormap.
%       Similar to clipping.
%
%  sgray(alpha)
%
%  IN  alpha: degree of BW color scale clustering(try 5)
%
%  Example: 
%
%  d = linear_events; imagesc(d); sgray(0.5);
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



i0=32;
i = 1:64;
t = (atan((i-i0)/alpha))';
s = t(64);
t = (t - min(t))*1./(max(t) -min(t));


m(1:64,2) = 1-t;
m(:,1) = 1-t;
m(:,3) =1-t;
colormap(m);
