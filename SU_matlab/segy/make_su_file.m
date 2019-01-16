function  [H] = make_su_file(filename,D,dtsec,offset,SCALCO);
%MAKE_SU_FILE: Write an su file to filename.
%
%  make_su_file(filename,D,dtsec,offset);
%
%  IN   filename: name of file where the data is saved
%       D:        data (traces are columns)
%       dtsec:    dt in secs
%       offset:   offset (integers)
% 
%  OUT  H:        header in structure 
%
%  This program only assigns very basic header words
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

[nt,nx]=size(D);
%cdp = 118:1950;

for k=1:nx,

 H(k)= make_empty_header;
 H(k).dt=floor(dtsec*1000*1000);
 H(k).d1=dtsec;
 H(k).ns=nt;
 H(k).trace =  D(:,k);
 H(k).tracl = k;
 H(k).offset = offset(k);
 H(k).scalco = SCALCO;
 %H(k).cdp = cdp(k);

end;

writesegy(filename,D,H);
