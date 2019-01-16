function [z] = zeros_wav(w);
%ZEROS_WAV: computes the zeros of a wavelet.
%
%  [z] = zeros(w);
%
%  IN   w: a wavelet;
%
%  OUT  z: complex zeros
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


[N,M] = size(w);
if N == 1; wr = fliplr(w); z = roots(wr); end;
if M == 1; wr = flipud(w); z = roots(wr); end;

return;
 

