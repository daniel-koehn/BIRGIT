function phase = calc_phase_error(t0,w,w1,dt,fmax,N,fphase)
%PLOT_SPECTRAL_ATTRIBUTES: Plot amplitude and phase of a wavelet.
%
%  plot_spectral_attributes(t0,w,dt,fmax,N,pos1,pos2);
%
%  IN     t0:    time in sec of the first sample of the
%                wavelet 
%         w:     input data or wavelet
%         dt:    sampling interval in secs
%         fmax:  max. freq. to display in Hz
%         N:     plot phase in the interval -N*180,N*180
%
%  Out    A figure in the current figure.
%         
%  Example: amplitude and phase of a 90 degree rotated wavelet
%
%         dt = 4./1000; fl=2; fh=90; c=90;
%         [w,tw] =rotated_wavelet(dt, fl, fh, c);
%         plot_spectral_attributes(min(tw),w,dt,125,1,221,222);
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
  FSize=20;
  nf = 4*2^nextpow2(length(w));
  n2 = nf/2+1;
  X = fft(w,nf); 
  X1 = fft(w1,nf); 
  fa = [0:1:nf-1]'/nf/dt;

% Tell the dft that the first sample is at t0 (to avoid unramping)

  Phase_Shift = exp(-1i*2*pi*fa*t0);
  X = X.*Phase_Shift;
  X1 = X1.*Phase_Shift;
  
  n2 = floor( fmax* (nf*dt)) +1;
  X = X(1:n2,1);
  X1 = X1(1:n2,1);
  
  f = [1:1:n2]'; f = (f-1)/nf/dt;
  Theta = angle(X.*conj(X1));
  
  Theta(1,1) = 0.;
  Theta = Theta*180/pi;
  nphase = round(fphase./f(2))+1;
  phase = Theta(nphase);
  
% Plot Phase

%     plot(f,Theta); xlabel('Frequency [Hz]'); ylabel('Phase [Deg]')
%     axis([0,fmax,-N*180,N*180]);
%     grid;
% 
%     set(get(gca,'title'),'FontSize',FSize);
%     set(get(gca,'title'),'FontWeight','bold');
%     set(get(gca,'Ylabel'),'FontSize',FSize);
%     set(get(gca,'Ylabel'),'FontWeight','bold');
%     set(get(gca,'Xlabel'),'FontSize',FSize);
%     set(get(gca,'Xlabel'),'FontWeight','bold');
%     set(gca,'FontSize',FSize);
%     set(gca,'FontWeight','bold');
%     set(gca,'Box','on');
%     set(gca,'Linewidth',1.0);            

