function [prim,m,tau,q] = pradon_demultiple(d,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cut);
%PRADON_DEMULTIPLE: Multiple removal using parabolic Radon Transform.
%
%  [prim,m,tau,q] = pradon_demultiple(d,dt,h,qmin,qmax,nq,flow,fhigh,mu,q_cut);
%
%  IN   d:     data   (d(nt,nh))
%       dt:    sampling interval in secs
%       h:     offset in meters  (h(nh))
%       qmin:  min residual moveout at far offset (secs)
%       qmax:  max residual moveout at far offset (secs)
%       nq:    number of samples of the residual moveout axis
%       flow:  min freq to process (Hz)
%       flow:  max freq to process (Hz)
%       mu:    trade-off parameter for LS solution used
%              to retrieve the Radon panel
%       q_cut: keep contributions from q_cut to qmax,
%              this is to estimate the multiples that 
%              are removed from the primaries
%
%  OUT  prim: primaries obtained by removing from the data the
%             multiples modelled  with the Radon transform
%       m:    panel with the Radon transform  (m(nt,nq))
%       tau:  vertical axis for the Radon panel (tau(nt))
%       q:    horizontal axis for the Radon panel (q(nq))
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination,
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%  Example: see radon_demo.m 
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


  [nt,nx] = size(d);
  dq = (qmax-qmin)/(nq-1);
  q = qmin+dq*[0:1:nq-1];

  N = 2;                % parabolic trasnsform

% Transform from t-x to tau-q

  [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,'ls');

  nq = length(q);
  iq_cut = floor((q_cut-qmin)/dq)+1;
  mc = m;

  mc(:,1:iq_cut) = 0;   % Keep multiples in the Radon panel

% Transform from tau-q to t-x

  [dm] = forward_radon_freq(mc,dt,h,q,N,flow,fhigh);

  prim = d-dm;

  tau = (0:1:nt-1)*dt;

  return;
