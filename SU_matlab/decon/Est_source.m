function [D_mod,D_s] = Est_source(nt,ntr,D_mod,D_data,D_spike,NORM_SOURCE,DAMP_SOURCE,out_wavelet,k,eps_w,max_offset_source,offset)
%   Est_source 
%   Estimate source wavelet by stabilized Wiener Deconvolution
%
%   Daniel Koehn
%   Kiel, the 3rd of January 2014

  for j=1:nt
            sumen(j)=0.0;
            sumed(j)=0.0;
            sumfn(j)=0.0;
            sumfd(j)=0.0;
        end
        
        for i=1:ntr
    
            % calculate FFT 
            D_model_fft = fft(D_mod(:,i));
            D_data_fft = fft(D_data(:,i));
   
            at = real(D_data_fft);
            bt = imag(D_data_fft);
   
            ct = real(D_model_fft);
            dt = imag(D_model_fft);
     
            a(:,i) = at;
            b(:,i) = bt;
   
            c(:,i) = ct;
            d(:,i) = dt;
      
        end
        
        % calculate nominator and denominator for Wiener filter
        for i=1:ntr    
            for j=1:nt
                    
                    if(abs(offset(i))<=max_offset_source)
                        
                        sumen(j) = sumen(j) + (a(j,i).*c(j,i)+b(j,i).*d(j,i));
                        sumed(j) = sumed(j) + (c(j,i).^2+d(j,i).^2);

                        sumfn(j) = sumfn(j) + (a(j,i).*d(j,i)-b(j,i).*c(j,i));
                        sumfd(j) = sumfd(j) + (c(j,i).^2+d(j,i).^2);
                        
                    end
                
            end
        end
        
        % calculate new source wavelet
        es = sumen./(sumed+eps_w.*max(abs(sumed)));
        fs = sumfn./(sumfd+eps_w.*max(abs(sumfd)));
        
        D_ss_fft = complex(es,fs);
        D_s_fft = D_ss_fft';
        
        D_spike_fft = fft(D_spike);
        %D_s  = real(ifft(D_s_fft.*D_spike_fft));
        D_s  = real(ifft(D_s_fft));
        
%         nth=round(nt./2.0);
%         
%         h=nth;
%         for i=1:nth
%             D_s1(i) = D_s(h);
%             h=h-1;
%         end
%         
%         h=nt;
%         for i=nth:nt
%             D_s1(i) = D_s(h);
%             h=h-1;
%         end
%         
%         D_s = flipud(D_s1');
        
        % calculate seismogram with new source wavelet
        D_s_fft = fft(D_s);
%         D_s_fft = D_s_fft.*D_spike_fft;
%         D_s = real(ifft(D_s_fft));
        
        for i=1:ntr
            D_model_fft = fft(D_mod(:,i));
            %D_data_fft = fft(D_data(:,i));
            D_sp = D_model_fft.*D_s_fft;
            D_mod(:,i) = real(ifft(D_sp));
        end
        
        if(NORM_SOURCE==1)
            D_s = D_s./max(abs(D_s));
        end
        
        if(NORM_SOURCE==2)
            
            D_s = D_s./abs(sum(D_s));
            %D_s = D_s./max(abs(D_s));
            
        end
        
        if(DAMP_SOURCE==1)
            
             % calculate discrete time samples
             nnorm0 = round(t0_source./DT);
             nnorm1 = round(t1_source./DT);
             
             % apply time damping at the beginning of the wavelet             
             D_s(1:nnorm0) = D_s(1:nnorm0).*exp(-adm_source.*(t(1:nnorm0)-t(nnorm0)).^2);
                  
             % apply time window after data traces
             if(nnorm1<=nt)
                D_s(nnorm1:nt) = D_s(nnorm1:nt).*exp(-adp_source.*(t(nnorm1:nt)-t(nnorm1)).^2);
             end
             
             clear nnorm0;
             clear nnorm1;
             
        end
     
        % output of estimated source wavelet
        out_wavelet1 = [out_wavelet,'_shot_',int2str(k),'.dat'];
        dlmwrite(out_wavelet1,D_s,'delimiter','\n','precision','%e');
        
        clear D_model_fft;
        clear D_data_fft;
   
        clear at;
        clear bt;
        clear ct;
        clear dt;

end

