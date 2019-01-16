function pimage(x,z,d);
%PIMAGE: High quality image for ppt presentations etc etc. 
%        Go to "file" "export setup" to export with proper labels, 
%        line width, size and other attributes.
%
%  IN  x: x-axis, x(nx)
%      z: z-axis, z(nz)
%      d: data, d(nz,nx)
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: U. Theune and M.D.Sacchi
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


 pcolor(x,z,d);
 shading flat;
 axis ij; 


 set(gca,'ydir','reverse','xaxislocation','top','yaxislocation','left','layer','top','linewidth',2,'tickDir','out','box','on')

 % Move the plot down to make space at the top

  pos=get(gca,'position');
  set(gca,'position',[pos(1) 0.05 pos(3) pos(4)])

 colormap(seismic(1));

 return;


