function [ D_data ] = norm_data(ntr, nt, D_data, NORM, NORM_RMS, NORMTR)
% Normalize modelled and field data
%
% Daniel Koehn
% Kiel, 11.11.2016

    % Data Normalization
    if(NORM==1) % Normalization to maximum value of each trace
        
        for i=1:ntr 
            
            if(NORM_RMS==0)
                max_data = max(abs(D_data(:,i)));
            end
            
            if(NORM_RMS==1)
                max_data = sqrt(sum(D_data(:,i).^2)./nt); 
            end
            
            % scale data amplitudes
            D_data(:,i) = D_data(:,i) ./ max_data;
         
        end
    end
    
    if(NORM==2) % Normalization to maximum value of each shot
            
            if(NORM_RMS==0)
                max_data = max(max(abs(D_data)));
            end
            
            if(NORM_RMS==1)
                
                for i=1:ntr
                    tmp(i) = sqrt(sum(D_data(:,i).^2)./nt); 
                end
                
                max_data = max(tmp);
                
            end
         
            % scale data amplitudes
            D_data = D_data ./ max_data;
            
    end
    
    if(NORM==3) % Normalization to maximum value of each shot
        
            if(NORM_RMS==0)
                max_data = max(abs(D_data(:,NORMTR)));
            end
            
            if(NORM_RMS==1)
                max_data = sqrt(sum(D_data(:,NORMTR).^2)./nt); 
            end
         
            % scale data amplitudes
            D_data = D_data./max_data;
         
    end

end

