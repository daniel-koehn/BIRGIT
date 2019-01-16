function [x] = laplace_mixture(N,arg);
%LAPLACE_MIXTURE: Deviates distributed according to a Laplace mixture.
%
% [x] = laplace_mixture(N,arg);
%
% IN   N:        Number of deviates
%      arg:      arg = [lambda1,lambda2,p] where 
%      lambda1:  Laplace paramater for the first distribution
%      lambda2:  Laplace paramater for the second distribution
%      p:        Mixing parameter (0,1)
%
% OUT  x(N,1): deviates
%
%
% Example with Walden and Hosken  parameters (Geophysical Prospecting, 1986):
%
% lambda1 = 0.007; lambda2 = 0.017; p = 0.24; 
% r = laplace_mixture(500, [lambda1, lambda2, p]);
% plot(r); title('Reflectivity');
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





 lambda1 = arg(1);
 lambda2 = arg(2);
 p = arg(3);

 x1 = laplace(N,lambda1);
 x2 = laplace(N,lambda2);

 x = zeros(N,1);
for k=1:N
  if  p>rand(1,1); x(k,1) = x1(k);
         else;
                   x(k,1) = x2(k);
   end
end

%expected_var = (p)*2*lambda1^2 + (1.-p)*2*lambda2^2
%var(x) 

return


