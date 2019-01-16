function plot_spectral_attributes_2(t0,data,tr1,tr2,dt,fmax,N,pos1,pos2);
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
%         pos1:  position of Ampl. in the figure (221)
%         pos2:  position of Phase in the figure (222);
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
  
  w = data(:,tr1);
  w1 = data(:,tr2);
  
  nf = 4*2^nextpow2(length(w));
  n2 = nf/2+1;
  X = fft(w,nf); 
  fa = [0:1:nf-1]'/nf/dt;

  X1 = fft(w1,nf);
  
% Tell the dft that the first sample is at t0 (to avoid unramping)

  Phase_Shift = exp(-i*2*pi*fa*t0);
  X = X.*Phase_Shift;
  n2 = floor( fmax* (nf*dt)) +1;
  X = X(1:n2,1);
  f = [1:1:n2]'; f = (f-1)/nf/dt;
  A = abs(X);
  A  = A/max(A);
  Theta = unwrap(angle(X));
  Theta(1,1) = 0.;
  Theta = Theta*180/pi;

  X1 = X1.*Phase_Shift;
  X1 = X1(1:n2,1);
  A1 = abs(X1);
  A1  = A1/max(A1);
  
% Plot Amplitue and Phase

  subplot(pos1); plot(f,A,'Linewidth',2.0);     
                 %xlabel('Frequency [Hz]'); 
                 ylabel('Amplitude');
                 title(['trace #',int2str(tr1)]);
                 axis([0,fmax,0,1.1])
                 grid;
       set(get(gca,'title'),'FontSize',FSize);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',FSize);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',FSize);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',FSize);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);

subplot(pos2); plot(f,A1,'Linewidth',2.0);     
                 xlabel('Frequency [Hz]'); 
                 ylabel('Amplitude');
                 title(['trace #',int2str(tr2)]);
                 axis([0,fmax,0,1.1])
                 grid;
       set(get(gca,'title'),'FontSize',FSize);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',FSize);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',FSize);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',FSize);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);