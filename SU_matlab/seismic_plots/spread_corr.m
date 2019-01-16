function [ D_data ] = spread_corr(ntr, nt, DT, fmax, w0, offset, t, D_data, D_mod, CORR, eps_corr, vph, offset_corr1, offset_corr2)
% Geometrical Spreading correction
% Implementation of different spreading corrections
%
% Daniel Koehn
% Kiel, 11.11.2016

  % apply spreading correction (Bleistein 1986) to field data
    % ---------------------------------------------------------
    if(CORR==1)
        for i=1:ntr
      
            tracei = D_data(:,i);
     
            % apply phase correction in frequency domain
            trace = spread_corr_phase(min(t),tracei,DT,fmax,w0,nt);
     
            D_data(:,i) = tracei(:) .* sqrt(t);

        end
    end
    
    % apply spreading correction, direct wave transformation (Forbriger et al 2014) to field data
    % -------------------------------------------------------------------------------------------
    if(CORR==2)                        
        
        for i=1:ntr
      
            tracei = abs(offset(i)) .* sqrt(2) .* D_data(:,i);            
            
            % apply phase correction in frequency domain
            trace = spread_corr_phase(min(t),tracei,DT,fmax,w0,nt);                    
     
            D_data(:,i) = tracei(:) .* sqrt(1.0./(t+eps_corr)); 
            
        end
    end
    
    % apply spreading correction sqrt(2 * offset * vph) to field data
    % ---------------------------------------------------------------
    if(CORR==3)
                        
        for i=1:ntr
                  
            tracei = D_data(:,i);
            
            % apply phase correction in frequency domain
            trace = spread_corr_phase(min(t),tracei,DT,fmax,w0,nt);                    
                 
            D_data(:,i) = tracei(:) .* sqrt(2.0 .* abs(offset(i)) .* vph);

        end
        
    end
    
    % apply geometrical spreading correction 
    % (hybrid transformation, Forbriger et al. 2014) to field data
    % ------------------------------------------------------------
    if(CORR==4)                
        
        taper_tmp = zeros(ntr,1);   
        
        % calculate taper function        
        D_OFFSET = abs(offset(2) - offset(1));
        ntaper = round(abs((offset_corr2 - offset_corr1)./D_OFFSET));
        ntaper2 = 2*ntaper;
        
        x = 1:ntaper2;
        % tmp =  sin(pi.*x./(ntaper2-1));
        tmp = 0.5 .* (1.0 - cos((2.0.*pi.*x)./(ntaper2)));
            
        for i=1:ntr
                                         
            % calculate taper function            
            if(i<=ntaper2)
                taper_tmp(i) = tmp(i);
            end
            
            if(i<=ntaper-1)
                taper_tmp(i) = tmp(ntaper);
            end
            
        end
        
        taper = taper_tmp(ntaper:ntr);
        
        clear taper_tmp;
        
        h = 1;
        for i=1:ntr
            
            tracei = D_data(:,i);             
            
            % apply phase correction in frequency domain
            tracei = spread_corr_phase(min(t),tracei,DT,fmax,w0,nt);                    
            
            % apply "single velocity transformation" for offset <= offset_corr1
            if(abs(offset(i)) <= offset_corr1)
                D_data(:,i) = tracei(:) .* sqrt(2.0 .* abs(offset(i)) .* vph);
            end
            
            % apply Cosine taper for offset_corr1 <= offset <= offset_corr2
            if((abs(offset(i)) > offset_corr1) && (abs(offset(i)) < offset_corr2))
                
                w1 = sqrt(2.0 .* abs(offset(i)) .* vph);
                w2 = abs(offset(i)) .* sqrt(2.0) .* sqrt(1.0./t);
                                            
                % apply taper function
                D_data(:,i) = tracei(:) .* (w1 .* taper(h) + w2 .* (1.0 - taper(h)));
                
                h = h + 1;
                clear tmp x;
                
            end

            % apply "direct wave transformation" for offset >= offset_corr2
            if(abs(offset(i)) >= offset_corr2)
                tracei = abs(offset(i)) .* sqrt(2.0) .* tracei;
                D_data(:,i) = tracei(:) .* sqrt(1.0./t);
            end                        
            
        end                
        
        %clear taper;
        
    end
    
    % apply geometrical spreading correction 
    % (polynomial fitting) to field data
    % --------------------------------------
    if(CORR==5)                
        
        % apply phase correction in frequency domain
        for i=1:ntr            
            tracei = D_data(:,i);                        
            trace = spread_corr_phase(min(t),tracei,DT,fmax,w0,nt);
            D_data(:,i) = trace;
        end
        
        % calculate RMS values of field data and modelled data
        for i=1:ntr
            AVO_data(i) = sqrt(sum(D_data(:,i).^2)./nt); 
            AVO_mod(i) = sqrt(sum(D_mod(:,i).^2)./nt); 
        end
        
        % assemble linear matrix equation
        % matrix A
        A = [ones(ntr,1) log(abs(offset))]
                
        % RHS vector
        b = log(AVO_mod./AVO_data)
        b = b';       
        
        % calculate solution vector
        x = (A'*A)\(A'*b)
        
        % apply polynomial correction to field data
        for i=1:ntr
            D_data(:,i) = exp(x(1)).*abs(offset(i)).^x(2).*D_data(:,i);  
        end
        
        
    end


end

