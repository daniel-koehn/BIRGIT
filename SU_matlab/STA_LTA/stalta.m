function varargout=stalta(sig,DT,BE,STA,LTA,TR,DTR,PEM,PET,PNL,ATL)
% [trigt,stav,ltav,ratio,tim1,tim2,tim3,trigs,dtrigs]=...
%   STALTA(sig,DT,BE,STA,LTA,TR,DTR,PEM,PET,PNL,ATL)
% 
% Triggering algorithm based on the ratio of the Short-Term
% Average absolute value to the Long-Term Average absolute value
% of a detrended signal. Everything is in standard units.
%
% INPUT:
%
% sig      Vector containing the signal
% DT       Sampling interval (s)
% BE       Beginning and end time of signal ([s s])
% STA      Short-term averaging window length (s)
% LTA      Long-term averaging window length (s)
% TR       Value of STA/LTA ratio that triggers
% DTR      Value of STA/LTA ratio that untriggers
% PEM      Time buffer added before triggering time (s)
% PET      Time-buffer added after detriggering time (s)
% PNL      Minimum window length of any triggered section (s)
% ATL      Time between trigger and detrigger that must be
%          exceeded in order for the triggered section to be reported (s)
%
% OUTPUT:
%
% trigt   Matrix with begin and end times of triggered sections (s)
% stav    Short-term average of absolute values of detrended signal
% ltav    Long-term average of absolute values of detrended signal
% ratio   Ratio of short-term to long-term average
% tim1    Time axis for 'ratio'
% tim2    Time axis for 'stav'
% tim3    Time axis for 'ltav'
% trigs   Vector with triggering points, in samples
% drtigs  Vector with detriggering points, in samples
%
%
% EXAMPLE:
%
% stalta('demo')
% stalta('demo',3)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/09/2009

%   % Default values
%   defval('DT',1)
% 
%   %defval('BE',[0 1]) % But this is no good, really
% % 
%   defval('STA',10)
%   defval('LTA',100)
%   defval('TR',2)
%   defval('DTR',1)
%   defval('PEM',100)
%   defval('PET',100)
%   defval('PNL',500)
%   defval('ATL',20)
  
  % Figure out how many samples the windows encompass
  STAsmp=ceil(STA/DT)
  LTAsmp=ceil(LTA/DT);
  PEMsmp=ceil(PEM/DT);
  PETsmp=ceil(PET/DT);
  NPTS=length(sig);

% Detrend the signal so DC value and trend don't play
 % [ %
 % Remove linear trend vom sig from data by finding the 
 % best fitting coeficients a,b, that describe
 % y=ax+b

trend=sig;

trend=trend(:);
m=length(trend);
x=1:m;

% Make Jacobian
A=[x(:) ones(m,1)];

% Invert it
% ab=inv(A'*A)*A'*trend;
[ab] = householder(A,trend);

a=ab(1);
b=ab(2);

% Calculate and remove linear trend
%sig=trend-a*x(:)-b;

%% disp('Detrending... may not be appropriate for everyone') 
  
% Calculate long-term and short-term average and their ratio
 % Moving average routine of signal abs(sig)  with window length 'XXXsmp'
 
y=abs(sig);
y = y(:);
[m n] = size(y);
if m < STAsmp
  disp('Window wider than sample');
  ma=0;
else
  stav=cumsum([sum(y(1:STAsmp));y(STAsmp+1:m)-y(1:m-STAsmp)])./STAsmp;
end

y=abs(sig);
y = y(:);
[m n] = size(y);
if m < LTAsmp
  disp('Window wider than sample');
  ma=0;
else
  ltav=cumsum([sum(y(1:LTAsmp));y(LTAsmp+1:m)-y(1:m-LTAsmp)])./LTAsmp;
end

  % If it's not possible, don't return anything
  if length(stav)==1 | length(ltav)==1
    [trigt,stav,ltav,ratio,tim1,tim2,tim3,trigs,dtrigs]=deal(NaN);
  else
    ratio=[repmat(stav(1),STAsmp-1,1) ; stav]...
	  ./[repmat(ltav(1),LTAsmp-1,1); ltav];
    % Simple solution to make last triggered point
    % automatically detriggered at the end of the section
    ratio(end)=0;
    
    %plot([repmat(stav(1),STAsmp-1,1) ; stav]);
    
    % Calculate timing axes for all three variables
    tim1=linspace(BE(1),BE(2),NPTS);
    tim2=linspace(BE(1)+STA-DT,BE(2),NPTS-STAsmp+1);
    tim3=linspace(BE(1)+LTA-DT,BE(2),NPTS-LTAsmp+1);

    % Trigger when ratio exceeds TR
    trigi=find(ratio>TR);
    if ~isempty(trigi) 
      cnt=1;
      trigs(cnt,1)=trigi(1);
      % Detrigger when ratio drops below DTR after the trigger
      dtrgi=find(ratio(trigi(1):end)<DTR)+trigi(1)-1;
      dtrigs(cnt,1)=dtrgi(1);
      
      % First trigger times
      trigt=[tim1(trigi(1))-PEM tim1(dtrgi(1))+PET];
      % After that (not in real time) need to work on the remainder
      trigi=find(ratio(dtrgi(1)+PETsmp:end)>TR)+dtrgi(1)+PETsmp-1;
      while ~isempty(trigi)
	cnt=cnt+1;
	trigs(cnt,1)=trigi(1);
	dtrgi=find(ratio(trigi(1):end)<DTR)+trigi(1)-1;
	dtrigs(cnt,1)=dtrgi(1);
	if ~isempty(dtrgi)
	  trigt(cnt,:)=[tim1(trigi(1))-PEM tim1(dtrgi(1))+PET];
	  trigi=find(ratio(dtrgi(1)+PETsmp:end)>TR)+dtrgi(1)+PETsmp-1;      
	else % Terminate loop
	  trigi=[];
	end
      end
    else
      trigt=[];
      trigs=[];
      dtrigs=[];
    end
    
    % Verify that the time between the real trigger and the end is 
    % at least PNL seconds long; if not, update
    for index=1:length(trigs)
      if (trigt(index,2)-[BE(1)+trigs(index)*DT])<PNL
	trigt(index,2)=[BE(1)+trigs(index)*DT+PNL];
      end
    end

    % Make sure none of the values exceed the data length or beginning
    trigt(trigt>BE(2))=BE(2);
    % But don't let it start at exactly zero, ever
    trigt(trigt<BE(1))=BE(1)+DT;

    % Make sure the actual trigger length is at least ATL
    kept=ceil([dtrigs-trigs]*DT)>=ATL;
    trigt=trigt(kept,:);
    trigs=trigs(kept,:);
    dtrigs=dtrigs(kept,:);
  end
  
  % Provide output
  vars={trigt,stav,ltav,ratio,tim1,tim2,tim3,trigs,dtrigs};
  varargout=vars(1:nargout);
