% Separate individual shots, definition of acquisition geometry
%
% Daniel Koehn
% Raisdorf, the 16th of January 2019
%
% SU matlab code partly from
% Copyright (C) 2008, Signal Analysis and Imaging Group
% http://www-geo.phys.ualberta.ca/saig/SeismicLab
% Author: M.D.Sacchi
%
% SegyMat codes from
% https://github.com/cultpenguin/segymat
% Author: Thomas Meier Hansen

clear all
close all

% Locate directories for SU matlab
Dir='../SU_matlab';
path(path, strcat(Dir,'/segy'));
path(path, strcat(Dir,'/SegyMAT'));
path(path, strcat(Dir,'/seismic_plots'));
path(path, strcat(Dir,'/bp_filter'));
addpath(Dir);

% Switch display types
% WIG = 1 (wiggle plot)
% WIG = 0 (image plot)
WIG = 0;

% clipping values
dtr = 2; % plot every dtr-th trace
clip=1.5e-1;
FSize=20; % font size
max_t = 1.0;

% define image plot clipping values
caxis_value2 = 2e1;
caxis_value1 = -caxis_value2;

% load colormap seismic
load 'seismic.map'

% define screen size
% Daniel:
screen_x = 1900;
screen_z = 1200;

% Eva:
% screen_x = 1900;
% screen_z = 1200;

% DELAY = 1 (apply time delay = on)
% DELAY = 0 (apply time delay = off)
DELAY=0;
delay=0.02;   % time delay for real data

% Normalization of field data
% NORM = 1 (normalization to maximum amplitude of each trace)
NORM=1;

% FILTER = 1 (trapez bandpass filter = on) 
% FILTER = 2 (Butterworth filter = on)
% FILTER = 0 (bandpass filter = off)
FILTER=2;

% corner frequencies
f1=0.0;
f2=5.0;
f3=25.0;
f4=50.0;

% Butterworth filter
order = 6;
fc = 20.0;

% basename of input SU file
infile=['raw_data/SH_shots_long_west/SH_long_west'];

% define acquisition geometry
minshot = 1;          % minimum shot no. to display
maxshot = 66;         % maximum shot no. to display
maxshot_true = 66;    % true maximum number of shots

sdatname = 'source/source_kleinneudorf_long_west.dat';     % specify name of source file 
rdatname = 'receiver/receiver_kleinneudorf_long_west.dat'; % specify name of receiver file 

% read receiver coordinates from receiver.dat
rec_data = load(rdatname);

% x,z-coordinates of the receiver
xrec = rec_data(:,1);
zrec = rec_data(:,2);

ntr=length(xrec);     % number of traces per shot

% clean up memory
clear rec_data;

% read source coordinates from source.dat
shot_data = importdata(sdatname, ' ', 1);

% x,z-coordinates of the source positions
xshot = shot_data.data(:,1);
zshot = shot_data.data(:,3);

% clean up memory
clear shot_data;

Fig = figure;
figure(Fig)

% loop over shots
for k=minshot:maxshot
    
    %if(k~=10 && k~=11 && k~=12 && k~=13 && k~=14 && k~=29)
    
    % load field data
    % ---------------
    infile1 = [infile,'_shot_',int2str(k),'.su'];
    [D_p,Htr,H] = ReadSu(infile1,'endian','b');
    
    % calculate time vector only for the 1st dataset
    if(k==minshot)
        
      DT=H(1).dt.*1e-6;    % get sample interval from SU-Header
      [nt,ntrg]=size(D_p); % number of time samples and traces
      t=1:nt;
      t=t.*DT;
      t=t';
        
    end
    
    if(DELAY==1)
       ndelay = round(delay./DT);
       zdelay = zeros(ntr,ndelay);
       D_p = [zdelay' ; D_p];

       tmp = D_p(1:nt,:);
       clear D_p;
       D_p = tmp;
       clear tmp;
    end
    
        % apply bandpass filter
    % ---------------------
    if(FILTER==1)  % simple trapez filter
        D_p =  bp_filter(D_p,DT,f1,f2,f3,f4);
    end
    
    if(FILTER==2)  % Butterworth filter
        D_p =  butter_filter(ntr,D_p,DT,fc,order,1);     
    end    
    
    % Normalize data to maximum-value of each shot
    if(NORM==2)
        for i=1:ntr
            D_p(:,i) = D_p(:,i)./max(max(abs(D_p(:,i))));
        end
    end
    
    offset = xrec;       % define offset range

    
    % plot data of individual shots
    % -----------------------------

    if(WIG==1)
        
        for i=1:dtr:ntr 
                 
            plot(offset(i)+(D_p(:,i)./(clip.*max(abs(D_p(:,i))))),t,'k-','Linewidth',4.0);
            hold on;
            
        end
    end

    % Display data as image plot 
    if(WIG==0)
         %colormap(gray);
         colormap(seismic);
         imagesc(offset,t,D_p);
         caxis([caxis_value1 caxis_value2]);
         hold on;
    end

           %set(gca,'DataAspectRatio',[1 250 1]);
           set(get(gca,'title'),'FontSize',FSize);
           set(get(gca,'title'),'FontWeight','bold');
           set(get(gca,'Ylabel'),'FontSize',FSize);
           set(get(gca,'Ylabel'),'FontWeight','bold');
           set(get(gca,'Xlabel'),'FontSize',FSize);
           set(get(gca,'Xlabel'),'FontWeight','bold');
           set(gca,'FontSize',FSize);
           set(gca,'FontWeight','bold');
           set(gca,'Box','on');
           set(gca,'Linewidth',1.0);
           axis ij
           axis([min(offset) max(offset) 0.0 max_t]);
           %axis equal
           %axis tight

           xlabel('xrec [m]');
           ylabel('time [s]');
           iter_text=['Kleinneudorf data shot no. ',int2str(k)];

           title(iter_text);

           set(gcf, 'PaperUnits', 'inches');
           set(gcf, 'PaperSize', [5 7]);

           set(Fig,'position',[0 0, screen_x screen_z])
           set(Fig,'PaperPositionMode','Auto') 
      
% pick first arrivals
[offsetTW1,tTW1]=ginput;
       
plot(offsetTW1,tTW1,'go',...
     'MarkerEdgeColor','g',...
     'MarkerFaceColor','g',...
     'MarkerSize',10);
       

      output=['pics/Kleinneudorf_psv_shot_',int2str(k)];
      %saveas(Fig,output,'psc2'); 
      saveas(Fig,output,'png');
       
% write picks of first arrivals to ASCII file   
imfile=['picks_TW2_',int2str(k),'.dat'];

fid = fopen(imfile,'w');
for i=1:length(offsetTW1)
      fprintf(fid,'%e \t %e \n',offsetTW1(i),tTW1(i));
end

fclose(fid);

    hold off;
    pause(0.5);
    clear D_p;
   
    % end
    
end
