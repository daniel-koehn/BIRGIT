function [f,o] = spiking(d,NF,mu);  
%SPIKING: Spiking deconvolution using Levinson's recursion.
%
%  [f,o] = spiking(d,NF,mu)
%
%  IN   d:  data (trace are columns)
%      NF: lenght of the spiking operator
%      mu: prewhitening in percentage  
%
%  OUT  f:  the filter
%       o:  the ouput or convolution of the data with 
%           the filter (adjusted to the length of the
%           input data and normalized).
%
%  Note: We assume a minimum phase wavelet, we also assume
%        that the reflectivity is a white process. The latter
%        allows us to estimate the autocorrelation of
%        the wavelet from the autocorrelation of the trace.
%
%  Reference: Robinson and Treitel, 1980, Geophysical Signal Analysis, Prentice Hall
%
%  Note: some clarity was lost in order to use Matlab function "levinson" 
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


 NF = NF - 1;
 [ns,ntraces] = size(d);
 dmax = max(max(abs(d)));

 R = xcorr(d(:,1),d(:,1),NF);       % Compute data autocorrelation for trace 1

 if ntraces>1;
  Ra = R;
   for k=2:ntraces;
    R = xcorr(d(:,k),d(:,k),NF);    % Compute data autocorrelation
    Ra = Ra + R;
   end;
  R = Ra/ntraces;
 end;

 Rs = R(:,1).*hamming(NF*2+1);
 r = Rs(NF+1:2*NF+1,1);
 r(1,1) = r(1,1)*(1 + mu/100.);     % Add pre-whitening for stability

 [f] = levinson(r,NF);              % Fast inversion of Toeplitz system 
                                    % via Levinson's recursive algorithm

 f = [f'];                          % I like column vectors

 if nargout == 2
  o = conv2(d,f);          
  o = o(1:ns,:);
  omax = max(max(abs(o)));
  o = o * dmax/omax;
 end

return

