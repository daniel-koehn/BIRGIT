
% EMD:  Emprical mode decomposition
%
% imf = emd(x,n)
%
% x   - input signal (must be a column vector)
% n   - number of intrinsic mode functions
%
% imf - Matrix of intrinsic mode functions (each as a row)
%
% See:  Huang et al, Royal Society Proceedings on Math, Physical, 
%       and Engineering Sciences, vol. 454, no. 1971, pp. 903-995, 
%       8 March 1998
%
% Author: Ivan Magrin-Chagnolleau  <ivan@ieee.org>
% 

function imf = emd(x,n);

c = x(:)'; % copy of the input signal (as a row vector)
N = length(x);

%-------------------------------------------------------------------------
% loop to decompose the input signal into n successive IMFs

imf = []; % Matrix which will contain the successive IMF, and the residue

for t=1:n % loop on successive IMFs
   
   %-------------------------------------------------------------------------
   % inner loop to find each imf
   
   h = c; % at the beginning of the sifting process, h is the signal
   SD = 1; % Standard deviation which will be used to stop the sifting process
   
   while SD > 0.3 % while the standard deviation is higher than 0.3 (typical value)
      
      % find local max/min points
      d = diff(h); % approximate derivative
      maxmin = []; % to store the optima (min and max without distinction so far)
      for i=2:N-2
         if d(i)==0                        % we are on a zero
            if sign(d(i-1))~=sign(d(i+1))  % it is a maximum
               maxmin = [maxmin, i];
            end
         elseif sign(d(i))~=sign(d(i+1))   % we are straddling a zero so
            maxmin = [maxmin, i+1];        % define zero as at i+1 (not i)
         end
      end
      
      if size(maxmin,2) < 2 % then it is the residue
         break
      end
      
      % divide maxmin into maxes and mins
      if maxmin(1)>maxmin(2)              % first one is a max not a min
         maxes = maxmin(1:2:length(maxmin));
         mins  = maxmin(2:2:length(maxmin));
      else                                % is the other way around
         maxes = maxmin(2:2:length(maxmin));
         mins  = maxmin(1:2:length(maxmin));
      end
      
      % make endpoints both maxes and mins
      maxes = [1 maxes N];
      mins  = [1 mins  N];
      
      
      %-------------------------------------------------------------------------
      % spline interpolate to get max and min envelopes; form imf
      maxenv = spline(maxes,h(maxes),1:N);
      minenv = spline(mins, h(mins),1:N);
      
      m = (maxenv + minenv)/2; % mean of max and min enveloppes
      prevh = h; % copy of the previous value of h before modifying it
      h = h - m; % substract mean to h
      
      % calculate standard deviation
      eps = 0.0000001; % to avoid zero values
      SD = sum ( ((prevh - h).^2) ./ (prevh.^2 + eps) );
      
   end
   
   imf = [imf; h]; % store the extracted IMF in the matrix imf
   % if size(maxmin,2)<2, then h is the residue
     
   % stop criterion of the algo. if we reach the end before n
   if size(maxmin,2) < 2
      break
   end
   
   c = c - h; % substract the extracted IMF from the signal
   
end

return
