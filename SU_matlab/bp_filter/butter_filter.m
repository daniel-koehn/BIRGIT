function  [o] =  butter_filter(ntr,d,dt,fc,order,sws)
%BP_Filter: Apply a band-pass filter to a group of traces.
%
% Butterworth-Filter
%
% Daniel Koehn
% Kiel, the 10th of January 2015
 
 % sampling rate
 fs = 1/dt;
 
 % low-pass filter
 if(sws==1)
    [B,A] = butter(order,fc/(fs/2),'low');
 end
 
 % high-pass filter
 if(sws==2)
    [B,A] = butter(order,fc/(fs/2),'high');
 end
 
 for k = 1:ntr
     tmp = d(:,k);
  o(:,k) = filter(B,A,tmp);
 end


 
