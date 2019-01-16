function r = bernoulli(N,arg);
%BERNOULLI:  Bernoulli deviates with parameter lambda and sigma.
%
% [r] = bernoulli(N,[lambda,sigma]);
%
% IN   N:    Number of deviates
%      arg = [lambda, sigma]; where 
%      lambda: Occurrence of a non-zero sample (0,1) 
%      sigma: standard error for non-zero samples 
%
% OUT  r(N,1): deviates
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


 lambda = arg(1);
 sigma = arg(2);

 r = zeros(N,1);

 for k=1:N
  if rand(1,1)>lambda; r(k,1)= sigma*randn(1,1); end
 end;




