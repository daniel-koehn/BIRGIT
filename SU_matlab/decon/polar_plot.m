function polar_plot(z);
%POLAR_PLOT: Plot the roots of a wavelet in polar coordinates.
%
%  polar_plot(z)
%
%  IN   z: zeros of the wavelet, they can be computed
%          with zeros_wav.m
%
%  NOTE: some zeros migth end up outside the plot (check axis)
%
%  Example
%
%    wavelet = ricker(20,0.004); z = zeros_wav(wavelet); polar_plot(z);
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

 
 x=cos(0:0.1:2*pi); y = sin(0:0.1:2*pi);
 
 II = find(abs(z)<1); z1=z(II);
 
 plot(x,y);hold on;plot(real(z),imag(z),'sk'); 
                  ;plot(real(z1),imag(z1),'*r'); 
axis equal;
 axis([-2,2,-2,2]);
return

