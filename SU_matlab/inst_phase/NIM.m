function [D_res] = NIM(D_mod,D_true,ntr,nt)
%   NIM 
%   Calculate non-integration method
%
%   Daniel Koehn
%   Kiel, the 26th of November 2015
       
    for i=1:ntr

        % calculate Q
        Qmod = 0.0;
        Qtrue = 0.0;
        for j=1:nt
            
            Qmod = Qmod + D_mod(j,i).^2; 
            Qtrue = Qtrue + D_true(j,i).^2; 
            Q_mod(j,i) = Qmod;
            Q_true(j,i) = Qtrue;
            
        end

        Q_mod(:,i) = Q_mod(:,i)./Qmod; 
        Q_true(:,i) = Q_true(:,i)./Qtrue;
        
        dQ(:,i) = Q_mod(:,i) - Q_true(:,i); 
        
        % calculate data residuals
        Q1 = 0.0;
        Q2 = 0.0;
        for j=1:nt
            
            Q1 = Q1 + dQ(j,i).*Q_mod(j,i); 
            Q2 = Q2 + Q_mod(j,i);
            
        end
        
        Q3 = 0.0;
        for j=nt:-1:1
           
            Q3 = Q3 + dQ(j,i);
            D_res(j,i) = -2.0.*D_mod(j,i).*(-Q3-Q1)./Q2;
                    
        end
        
    end
    
end



