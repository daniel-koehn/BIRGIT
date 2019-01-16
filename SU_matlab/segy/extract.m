function [value] = extract(filename,hw)
%EXTRACT: Extract a header word from a SU-SEGY file. 
%
%   example:    cdp = extract('data','cdp') will provide
%               the cdp numbers of each trace in 'data' 
%               See documentation to see the rest of the header
%               words
%
%   example:    of = extract('data','offset') will extract offset
%               of each trace
%
%   example:    D = extract('data','trace') will extract the traces
%               into a matrix called D
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
~

   FID = fopen(filename,'r','l');      % l is for little endian

   segy=segy_struct;                   % load the definitions of
                                       % the header words

   count=count_struct;                 % load the position of   
                                       % each word in the header (in bytes)

   status = fseek(FID,count.ns,'bof'); % go to begging of file

   ns = fread(FID,1,segy.ns);          % read ns from first trace

   if strcmp(hw,'ns') == 0;            % get only ns and return
    else
     value = ns;
      return
       end

   total = 60+ns;                      % total num of 4bytes words 
                                       % in the header

   if strcmp(hw,'trace') == 0;              
    elements = 1;
     word=zeros(1,1);
      else
      elements = ns;
       word=zeros(ns,1);
        end

   hc=eval(strcat('count.',hw));
   hp=eval(strcat('segy.',hw));


   max_traces = 999999;

   for j=1:max_traces;

     position = total*(j-1)*4+hc;
      status = fseek(FID,position,'bof'); 

       if status == 0                  % check status and read 
        word = fread(FID,elements,hp);
        j = j + 1  
        value(:,j-1)  = word(:,1);
      else
       return                          % return if end of file reached 
      end 
   end 

 [message,errnum] = ferror(FID)



   

