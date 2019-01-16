function [houseresult] = householder(H, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Linear Least Squares with Householder method
%
%   householder benutzt die Householder Transformation zur Loesung des
%   
%   linearen, ueberbestimmten Gleichungssystems.                      
%   x ist die Loesung des Minimierungsproblems mat * x - b im Sinne   
%   der euklidschen Norm. Diese Loesung des Minimierungproblems muss  
%   nicht notwendigerweise auch eine Loesung des Gleichungssystems    
%   sein (Pruefung !)   
%
%    Eingabeparameter:                                                
%   ================                                                                           
%      H             Matrix des Gleichungssystems. H(m, n)                      
%                               
%      x             Rechte Seite des Gleichungssystems. x(m)                  
%                                                                                                                                                  
%    Ausgabeparameter:                                                
%   ================                                                 
%      houseresult   Loesungsvektor des Systems.    houseresult(m)                                      
%                                         
%                                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOUSEHOLDER Method   ==>    Test mit matlab Befehl x = H\z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                     
m =size(H,1);
n =size(H,2);  
mat = H;
b=z(:,1);

  for i=1:n %(i = 0; i < n; i++)       %  Householder Transformation   
    r = 0.0;
    for k=i:(m) %(k = i; k < m; k++)
      r = r + mat(k, i) * mat(k, i);
    end
    
%     if (r == 0.0)                    %  Matrix hat nicht Hoechstrang 
%       return (2);
%     end
        
    if mat(i, i) >= 0.0
      alpha = sqrt(r);
    else
      alpha = - sqrt(r);                                   
    end                           
    ak = 1 / (r + alpha * mat(i, i));
    mat(i, i) = mat(i, i) + alpha;

    d(i) = - alpha;

    maxnorm = 0.0;
    for k = (i + 1):(n) %(k = i + 1; k < n; k++)
 
      norm = 0.0;
      f = 0.0;
      for j = i:(m) %(j = i; j < m; j++)
        tmp = mat(j, k);
        f = f + tmp * mat(j, i);
        norm = norm + tmp * tmp;
      end

      if norm > maxnorm                 
        maxnorm = norm;
      end 
      
      f = f*ak;
      for j = i:(j) %(j = i; j < m; j++)
        mat(j, k) = mat(j, k) - f * mat(j, i);
      end
    end
%     if (abs(alpha) < eps * sqrt(maxnorm))    % Loesbar ?          
%     {
%       return (3);
%     }

    f=0;
    for j = i:(m)  %(f = 0, j = i; j < m; j++) % Rechte Seite transformieren  
      f = f + b(j) * mat(j, i);
    end
    
    f = f * ak;
    for j = i:(m)  %(j = i; j < m; j++)
      b(j) = b(j) - f * mat(j, i);                
    end
    
  end %for i

  for i = fliplr(1:n)  %(i = n - 1; i >= 0; i--) % Rueckwaertselimination       

    sum = 0;
    for k = (i + 1):(n) % (k = i + 1; k < n; k++)
      sum = sum + mat(i, k) * b(k);
    end
    
  b(i) = (b(i) - sum) / d(i);
  end
 
  
houseresult = b(1:n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%