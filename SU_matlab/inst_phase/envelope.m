function [D_res] = envelope(D_mod,D_true,ntr,nt)
%   inst_phase 
%   Calculate instantaneous phase data residuals
%
%   Daniel Koehn
%   Kiel, the 4th of November 2015

    % calculate Hilbert transform        
    for i=1:ntr

        % Hilbert transform
        D_mod_ana = hilbert(D_mod(:,i));
        D_true_ana = hilbert(D_true(:,i));

        E_mod = sqrt((real(D_mod_ana)).^2+(imag(D_mod_ana)).^2);
        E_true = sqrt((real(D_true_ana)).^2+(imag(D_true_ana)).^2);
        
        % Envelope of modelled data
        for j=1:nt
           
            D_res(j,i) = 0.0;
            
            if(abs(D_true(j,i))>0.0)
              D_res(j,i) = E_mod(j) - E_true(j);
            end
        
        end
        
    end

end

