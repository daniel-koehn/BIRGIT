function  [w,tw] =  trapezoidal_wavelet(dt,f1,f2,f3,f4,c);
%TRAPEZOIDAL_WAVELET: Computes a band-pass wavelet with phase rotation c.
% 
%  [w,tw] = make_a_trapezoidal_wavelet(dt,f1,f2,f3,f4,c);
%
%  IN     dt:   sampling interval in sec
%         f1:   freq. in Hz
%         f2:   freq. in Hz
%         f3:   freq. in Hz
%         f4:   freq. in Hz
%         c:    rotation in degs
%
%   ^
%   |     ___________
%   |    /           \   Amplitude spectrum
%   |   /             \
%   |  /               \
%   |------------------------>
%      f1 f2        f3 f4
%
%  OUT    w:    wavelet (column)
%         tw:   time axis in secs
%
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



 fc = (f3-f2)/2;
 L = floor(1.5/(fc*dt));
 nt = 2*L+1;
 k = nextpow2(nt);
 nf = 4*(2^k);

 i1 = floor(nf*f1*dt)+1;
 i2 = floor(nf*f2*dt)+1;
 i3 = floor(nf*f3*dt)+1;
 i4 = floor(nf*f4*dt)+1;

 up =  (1:1:(i2-i1))/(i2-i1);
 down = (i4-i3:-1:1)/(i4-i3);
 aux = [zeros(1,i1), up, ones(1,i3-i2), down, zeros(1,nf/2+1-i4) ];
 aux2 = fliplr(aux(1,2:nf/2));

 F = ([aux,aux2]');
 Phase = (pi/180.)*[0.,-c*ones(1,nf/2-1),0.,c*ones(1,nf/2-1)];
 Transfer = F.*exp(-i*Phase');
 temp = fftshift((ifft(Transfer)));
 temp =real(temp);

 
 w = temp(nf/2+1-L:nf/2+1+L,1);
 nw = length(w);
 w = w.*hamming(nw);

 if nargout>1;
  tw = (-L:1:L)*dt;
 end;


 

