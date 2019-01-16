function [ Data ] = agc( Data, window, DT )
% MATLAB 5.0 script for applying AGC to CSG traces in a CSG. 
% For each trace, the energy within a moving window of length
% np is computed and the sample in the center part of this window
% is divided by this energy.
%
% author: Gerhard Schuster
%
% data(x,t)   - input - nx x nt CSG matrix 
% np          - input - length of AGC window
% out1        - output- traces with AGC applied
%

  for it=1:size(Data,2)
      
%     disp([mfilename,' AGC trace : ',num2str(it)])
        
    % Trace Data
    TraceData=Data(:,it);
    
    % Window Length
    nsw=round(window./DT);
    nshalf=floor(nsw/2);
    startsample=ceil(nsw/2);
    endsample=length(TraceData)-floor(nsw);
    
    for is=startsample:1:endsample;
      
      range=[is-nshalf+1:1:is+nshalf];
      
      gain=mean(abs(TraceData(range)));
      if gain~=0
	
        TraceData(is)=TraceData(is)./mean(abs(TraceData(range)));	

        % APPLY TO TOP
        if is==startsample
          for i=[1:startsample-1];
            TraceData(i)=TraceData(i)./mean(abs(TraceData(range)));
          end
        end

        if is==endsample
          for i=[endsample+1:1:length(TraceData)];
            TraceData(i)=TraceData(i)./gain;
          end
        end
	
      end
      
    end
    
    Data(:,it)=TraceData;
    
  end
  
end

