function [D_data] = stat_corr(rec_stat,DT,ntr,nt,D_data,mode)
% Static correction
% apply static corrections to seismic data
% 
% Daniel Koehn
% Kiel, the 5th of October 2016

    ndelay = round(rec_stat./DT);
    tmp = zeros(nt,ntr);

    % apply receiver static correction
    if(mode==1)

       for i=1:ntr

           % ndelay(i)      

           for j=1:nt
               if((j+ndelay(i)>=1)&&(j+ndelay(i)<=nt))
                   tmp(j,i) = D_data(j+ndelay(i),i);
               end
           end

       end

    end

    % apply source static correction
    if(mode==2)
        
       for i=1:ntr

           % ndelay      

           for j=1:nt
               if((j+ndelay>=1)&&(j+ndelay<=nt))
                   tmp(j,i) = D_data(j+ndelay,i);
               end
           end

       end

    end
    
    D_data = tmp;
    clear tmp;
    clear ndelay;

end

