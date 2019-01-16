function [x] = laplace(N,lambda);
%LAPLACE: Laplace deviates with parameter lambda.
%         Distribution  f(x) = 1/2/lambda * e{-|x|/lambda}}
%
% [x] = laplace(N,lambda);
%
% IN   N:      Number of deviates
%      lambda: Laplace parameter
%
% OUT  x(N,1): deviates
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
%     



 y = rand(N,1);

 for k=1:N;
   if y(k)<=0.5; x(k)= lambda*log(2*y(k)     );    else
                 x(k)=-lambda*log(2*(1-y(k)));
  end;
 end;

 return;
