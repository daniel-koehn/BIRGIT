function [D_res] = inst_phase(D_mod,D_true,ntr)
%   inst_phase 
%   Calculate instantaneous phase data residuals
%
%   Daniel Koehn
%   Kiel, the 4th of November 2015

    eps_ip = 1e2;
    eps_E = eps_ip;

    % calculate Hilbert transform        
    for i=1:ntr

        % analytical signals
        D_mod_ana(:,i) = hilbert(D_mod(:,i));
        D_true_ana(:,i) = hilbert(D_true(:,i));

        % Envelope of modelled data
        E(:,i) = sqrt((real(D_mod_ana(:,i))).^2+(imag(D_mod_ana(:,i))).^2);

    end

    D_mod_ana_max = max(max(real(D_mod_ana)));
    D_true_ana_max = max(max(real(D_true_ana)));
    E2_max = max(max(E.^2));

    for i=1:ntr

        % instantaneous phase
        phi_mod = atan(imag(D_mod_ana(:,i))./(real(D_mod_ana(:,i))+eps_ip.*D_mod_ana_max));
        phi_true = atan(imag(D_true_ana(:,i))./(real(D_true_ana(:,i))+eps_ip.*D_true_ana_max));

        % Hilbert of (phi_true-phi_mod).*D_mod_ana./E.^2
        E_hilb = imag(hilbert((phi_true-phi_mod).*D_mod_ana(:,i)./(E(:,i).^2 + eps_E.*E2_max)));

        % Hilbert transform of analytical signals
        D_mod_ana_hilb = imag(hilbert(D_mod_ana(:,i)));

        % calculate data residuals   
        D_res(:,i) = (phi_true-phi_mod).*D_mod_ana_hilb./(E(:,i).^2 + eps_E.*E2_max) + E_hilb;

    end

end

