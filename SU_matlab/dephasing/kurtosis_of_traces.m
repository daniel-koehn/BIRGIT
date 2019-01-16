function [K] = kurtosis_of_traces(x)
%KURTOSIS: Kurtosis of a one or more time series.
%
%  [K] = kurtosis_of_traces(x);
%
%  IN   x:      data (vector or matrix)
%
%  OUT  K:      Kurtosis
%
%
%  Kurtosis is defined as K = E(x^4)/ (E(x^2))^2;
%
%    K = 3 for a Gaussian series 
%
%  There is a different definiton k' = K-3
%
%    k' = 0 for a Gaussian series
%
%  k' is called the Kurtosis Excess. This matlab 
%  function computes K not k'.
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


 N = ndims(x);

[n1,n2]=size(x);
 Ns = n1*n2;

 if N==1; 
  sum1 = sum(x.^4)/Ns;
  sum2 = (sum(x.^2)/Ns)^2;
 end;
 
 if N==2; 
  sum1 = sum(sum(x.^4)/Ns);
  sum2 = (sum(sum(x.^2))/Ns)^2;
 end;

 K = sum1/sum2;

 return;
