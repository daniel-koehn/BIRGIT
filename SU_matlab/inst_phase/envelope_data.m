function [D_true] = envelope_data(D_true,ntr,nt)
%   inst_phase 
%   Calculate envelope data
%
%   Daniel Koehn
%   Kiel, the 4th of November 2015

    % calculate Hilbert transform        
    for i=1:ntr

        % Hilbert transform
        D_true_ana = hilbert(D_true(:,i));
        E_true = sqrt((real(D_true_ana)).^2+(imag(D_true_ana)).^2);
        
        % Envelope of modelled data
        for j=1:nt           
           D_true(j,i) = E_true(j);
        end
        
    end

end

