function  plot_wb(t,w,a);
%PLOT_WB: Horizontal plotting with bias.
%
%  IN    t: time axis
%        w: traces or wavelets in columns.
%        a: overlab in percentage
%
%
%  Example:
%            
%      [w1,tw] = rotated_wavelet(4./1000,4.,50.,0);
%      [w2,tw] = rotated_wavelet(4./1000,4.,50.,45);
%      [w3,tw] = rotated_wavelet(4./1000,4.,50.,90);
%      W =[w1,w2,w3];
%      subplot(221); plot_wb(tw,W,  0);
%      subplot(222); plot_wb(tw,W,50);
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

 w = w/max(max(abs(w)));
 [nt,nx] = size(w);
 x = 1:1:nx;
 plot(t,((a/100+1)*w)+ones(size(w))*diag(x),'b','LineWidth',1.);axis tight;
 return;
