function [f,Q,U,As] = Q_spec_div(d,ds,dt,fmax,offset)
%FK_SPECTRA: FK spectrum of a seismic gather.
%
%  [S,k,f] = fk_spectra(d,dt,dx,L);
%
%  IN   d:      data (traces in columns) 
%       dt:     time interval
%       dx:     spatial increment between traces 
%       L:      apply spectral smoothing using a separable
%               2D Hamming window of LxL samples
%
%  OUT  S:      FK spectrum
%       f:      freq axis in Hz
%       k:      wave-number spectrum in cylces/m (if dx is in meters)
%
%  Note: when plotting spectra (S)  use log(S) or S.^alpha (alpha=0.1-0.3) to
%        increase the visibility of small events 
%
%  Example: 
%
%    [d,h,t] = linear_events; dt = t(2)-t(1); dx = h(2)-h(1); 
%    [S,k,f] = fk_spectra(d,dt,dx,6); imagesc(k,f,S);
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
%     

  [nt,nx] = size(d);

  nf = 4*2^nextpow2(nt);
  n2 = nf/2+1; 
  fa = [0:1:nf-1]'/nf/dt;

  cphase = 64.0;
  
 for i=1:nx
     
      tmp = d(:,i);
      X = fft(tmp,nf);

      Phase_Shift = exp(-i*2*pi*fa*dt);
      X = X.*Phase_Shift;
      n2 = floor( fmax* (nf*dt)) +1;
      X = X(1:n2,1);
      f = [1:1:n2]'; f = (f-1)/nf/dt;
      A = abs(X);
      U(:,i) = A;
     
 end
 
 % fft of source wavelet
 Xs = fft(ds,nf);
 Xs = Xs.*Phase_Shift;
 Xs = Xs(1:n2,1);
 As = abs(Xs);
 
 clear tmp A;
 
 % Qs-estimation by spectral division
 % divide amplitude spectra of trace i by trace j 
 for j=6:6
     h=1;
     for i=21:21
         
         Q(:,h) = -2.0.*pi.*f.*(offset(i)-offset(j))./(2.0.*cphase.*log(offset(i).*U(:,i)./(offset(j).*U(:,j))));
         h=h+1;
         
     end
     
 end
 

 return;
