function [D,H] = readsegy(filename,hw,min,max)
%READSEGY: Read SU-SEGY data.
%          The data and headers can be  extracted  in 
%          a range given by hw (header words).
%   
%  [D,H] = readsegy(filename,hw,min,max) returns the data D and
%
%  IN   filename: name of segy
%       hw:       header word to limit number of traces to read
%       min, max: min and maximum value of hw
%
%  OUT  D:        the data (in a matrix) 
%       H:        the header in a structure
%
%
%  Examples:    
%
%    [D,H] = readsegy('data_with_noise.su');
%    [nt,nh] = size(D);
%    offset = [H.offset];
%    dtsec = H(1).dt/1000/1000;
%    imagesc(offset,[0:1:nt-1]*dtsec,D);
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

 
  FID = fopen(filename,'r','l');       % l for little endian

   segy=segy_struct;                   % load the definitions of 
                                       % the header words

   count=count_struct;                 % load the position of 
                                       % each word in the header (in bytes)

   status = fseek(FID,count.ns,'bof'); % go to the beggining of file

   ns = fread(FID,1,segy.ns);          % read ns from first trace
                                       % ns is the number of traces per trace

   total = 60+ns;                      % total nuber of 4bytes words


   max_traces=9999999;                 % maximum number of traces (will
                                       % stop before). The variable status
                                       % will make the code stop when
                                       % the end of file is reached  

 if nargin~=1;
   hc=eval(strcat('count.',hw));       % assigned the header word required  
   hp=eval(strcat('segy.',hw));        % to extract the traces.
 j = 1;                                % counter 
 for k =1:max_traces
   position = total*(k-1)*4+hc;        % where in bytes is the header word
   status = fseek(FID,position,'bof'); 
  if status == 0                       % stop when end of file is reached
    w = fread(FID,1,hp);
     if  w>=min                        % pick traces with a given range
      if w<=max                        % of the desired header word
       position = total*(k-1)*4+count.trace;
       status = fseek(FID,position,'bof'); 
       trace = fread(FID,ns,segy.trace);
       j = j + 1;  
       D(:,j-1)  = trace(:);           % load traces into D
       H(j-1)  = header(FID,ns,k);     % load each header in a structure H
%       disp(j-1)
      end
     end
   else
   sprintf(' Number of traces =  %d',j-1)
  return
 end
 end 


 else 

% when no hw and bounds are given, reads evrything

 for k =1:max_traces
       position = total*(k-1)*4+count.trace;
       status = fseek(FID,position,'bof'); 
       if status == 0                 
       trace = fread(FID,ns,segy.trace);
       D(:,k)  = trace(:);           % load traces into D
       H(k)  = header(FID,ns,k);     % load each header in a structure H
%      disp(k)
   else
   sprintf(' Number of traces =  %d',k-1)
  return
  end
 end 
Ntraces=k-1;
 end 



 [message,errnum] = ferror(FID)
 fclose(FID);

