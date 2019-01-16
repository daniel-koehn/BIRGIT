function [refl,dp] = sparse_decon(d,w,mu,iter_max);
%SPARSE_DECON: Sparse-spike deconvolution using a l1 norm regularization.
%              Find refl by minizing ||w*refl - d||^2 + mu l1(refl)
%
% [refl,dp] sparse_decon(d,w,mu,iter_max);
%
% IN   d(nt,1):  trace
%      w(nw,1):  wavelet (I delay output of decon assuming
%                a zero phase wavelet is used)
%      mu:       regularization parameter 
%      iter_max: maximum number of IRLS iterations
%
% OUT  refl:     broad-band reflectivity
%      dp:       predicted trace (convolution of refl with wavelet)
%
%
%  Note: l1 regularization with Iterative Reweigthed Least Squares (IRLS)
%
%  Reference: Sacchi, M.D., 1997, Re-weighting strategies in seismic 
%             deconvolution: Geophysical Journal International, 129, 651-656.
%
%  Example: see sparse_decon_demo.m
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

  nw = length(w);
  [nt,ntraces] = size(d);

  d = [d;zeros(nw-1,ntraces)];
  W = convmtx(w,nt);
  R = W'*W;                           % Autocorrelation of wavelet
  R0 = trace(R);                      % Trace of Autocorrelation

  Q = mu*R0*eye(nt);                  % Initial Matrix of Weights
  sc = 0.0001;

 for itrace = 1:ntraces;

  s = d(:,itrace);
  r = zeros(nt,1);                    % Initial solution

  for k=1:iter_max;                   % Start Iterative algorithm

    g = W'*s;
    Matrix = R + Q;
    r = Matrix\g;
    Q = mu*diag(1./(abs(r)+sc));

  end;

    refl(:,itrace) = r;

 end;

 n2 = floor(nw/2)
 temp = refl;
 dp = conv2(refl,w);
 dp = dp(1:nt,:);
 refl = [zeros(n2,ntraces);temp(1:nt-n2,:)];


 return;
