function [P,f,w,tw]  = smooth_spectrum(d,dt,L,io);
%SMOOTH_SPECTRUM: Power spectrum estimate by smoothing the periodogram.
%                 For more than one trace provides the average spectrum
%                 followed by smoothing.
%
%  [P,f,w,tw]  = smooth_spectrum(d,dt,L,io);
%
%  IN     d:  data (traces are columns)
%         dt: sampling interval in secs
%         L:  Lenght of the freq. smoothing operator
%             L=0  means no snoothing 
%         io: 'db' in db, 'li' for linear scale 
%
%  OUT    P:  normalized smoothed power spectrum in linear or dB scale.
%         f:  frequency axis 
%         w:  wavelet (zero phase with Power spectrum P)
%         tw:  time axis for wavelet
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


 [nt] = size(d,1);
 wind = hamming(2*L+1);

 nf = max(2*2^nextpow2(nt),2048);

 f = 1:1:nf/2+1;
 f = (f-1)/dt/nf;      % Freq. axis in Hz 

 D = fft(d,nf,1);

 D = abs(D).^2;


 ND = ndims(d);
 

  if ND==2; 
  [nt,nx]=size(d) ;  D =  sum(D,2)/nx; 
  end;

  if ND==3; 
  [nt,nx,ny]=size(d);  
  nx
   ny
D =  sum(squeeze(sum(D,3)),2)/(nx*ny); 

 
  end;

 D = conv(D,wind);     % Smooth 
 N = length(D);
 D = D(L+1:N-L);
 A = sqrt(D);

 D = D(1:nf/2+1,1);
 D = D/max(D);

if(io=='db'); 
 P = 10*log10(D);      % Power spectrum in dB
 I = find(P<-20); P(I)=-20;
else
P = D;
end

if (nargout>2); 

% Make a wavelet    
% For the length I am assuming a wavelet with central
% freq. 30Hz

  f0 = 30;

  L = 3/(f0*dt);
  w = real(fftshift( ifft(A)));
  w = w(nf/2+1-L:nf/2+1+L,1);
  w = w.*hamming(2*L+1);
  w = w/max(abs(w));
  tw = [-L:1:L]*dt;

end;

