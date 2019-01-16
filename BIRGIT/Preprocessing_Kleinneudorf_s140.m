% Separate individual shots, definition of acquisition geometry
%
% Daniel Koehn
% with contributions by ...
% Eva Dokter
% Denise De Nil
% Kiel, the 16th of January 2019
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

% Locate directories for SU matlab/SegyMat
Dir='../SU_matlab';
path(path, strcat(Dir,'/segy'));
path(path, strcat(Dir,'/spectra'));
path(path, strcat(Dir,'/SegyMAT'));
path(path, strcat(Dir,'/seismic_plots'));
path(path, strcat(Dir,'/bp_filter'));
path(path, strcat(Dir,'/fk_filter'));
path(path, strcat(Dir,'/STA_LTA'));
path(path, strcat(Dir,'/radon_transforms'));
path(path, strcat(Dir,'/phase_vel_freq'));
path(path, strcat(Dir,'/decon'));
path(path, strcat(Dir,'/inst_phase'));
addpath(Dir);

WRITE_ACQ = 0;

% WRITE_DATA = 1 (write frequency data and pick files = on)
% WRITE_DATA = 2 (write time domain data and pick files = on)
% WRITE_DATA = 0 (write data and pick files = off)
WRITE_DATA = 0;

% define frequencies
nf = 240;
fmin = 10.0;
fmax = 20.0;
df = (fmax - fmin) / nf; 
f = fmin:df:fmax;

% define parameters for fk-filter
FK_FILT = 0;
fmaxfk = 50.0;
Lham = 20;
vmin = 50.0;
vmax = 100.0;

% damp surface waves in time domain
DAMP_SW = 0;
maxoffset_sw = 4000.0;
gamma_sw = 10.0;
%vmin_sw = 0.0;
%vmax_sw = 200.0;

vmin_sw1 = 0.0;
vmax_sw1 = 40.0;

% calculate envelope of data
ENVELOPE = 0;

% Switch display types
% WIG = 1 (wiggle plot - mod+field data / data residuals)
% WIG = 0 (image plot - field data)
% WIG = 2 (image plot - field data / model data)
% WIG = 3 (AVO plot - field data / model data)
% WIG = 4 (wiggle plot - field data only)
% WIG = 5 (plot phase error difference)
% WIG = 6 (phase velocity - frequency spectra)
% WIG = 7 (Q-value estimation by spectral division)
% WIG = 8 (amplitude and phase spectra residuals)
% WIG = 9 (fk-spectra)
% WIG = 10 (plot acquisition geometry)
% WIG = 11 (stack data)
% WIG = 12 (wiggle plot - mod+field data)
WIG = 12;

% phase velocities used for plotting phase velocity-frequency spectra
cphase = 20.0:0.5:500.0;

% parameters for phase error calculations
FWI_phase_error=0;
fphase=10.0;      % frequency for phase error calculation

% parameters for WIG=8
tr_freq = 10;   % calculate amplitude spectra for trace tr_freq 

% Define data residuals
% MISFIT = 2 (L2-Norm)
% MISFIT = 5 (Gobal correlation norm)
MISFIT = 2;

% Exponential gain of data
GAINDATA = 0;
again = 3e0;

% apply static correction
STAT_REC = 0;
STAT_SRC = 0;

% clipping values
dtr = 1; % plot every dtr-th trace
dresamp_nt = 2;  % resample wiggle plot to reduce image size
clip=1e0;
FSize=20; % font size
Tmin=0.1;
Tmax=-0.8;

% define image plot clipping values
caxis_value2 = 2e-3; % NORM=2
% caxis_value2 = 1e-6; % NORM=1
caxis_value1 = -caxis_value2;

% load colormap seismic
load 'seismic.map'

% define screen size
% full dataset
screen_x = 1900;
screen_z = 1200;

% DELAY = 1 (apply time delay = on)
% DELAY = 0 (apply time delay = off)
DELAY=1;
delay=0.1;         % time delay for real data
gamma_delay = 1e0;   % time delay damping

% name of output SU data file
out_su_file=['output_data/DENISE_kleinneudorf_SH'];
SCALCO = -1000;

% Normalization of field data
% 
% NORM = 1 (normalization to maximum amplitude of each trace)
% NORM = 2 (normalization to maximum amplitude of each shot)
% NORM = 3 (normalization to maximum amplitude of each trace with respect to each dataset)
NORM=2;

% TWIN = 1 (time window = on)
% TWIN = 0 (time window = off)
% TWIN = 2 (read picked times)
% TWIN = 3 (read picked times + use constant twin+ at constant time)
% TWIN = 5 (STA/LTA picker)
TWIN=0;

% Parameters for TWIN=1 (linear window)
twinc=0.0;     % center of time window
twinp=0.07;  % twinc + twinp
twinm=0.0;    % twinc - twinm
vnmo=1500;     % velocity of the Rayleigh wave

% General parameters for TWIN>0
% define spline fit
delay_tw = 0.0; % time delay of twin- und twin+

adm=1e7;          % damping parameter before first arrival
adp=1e7;          % damping parameter after end of time window

twinp_const=1.6;  % if TWIN==3 use constant twin+ at constant time 

% TDAMP = 1 (zero-time damping = on)
% TDAMP = 0 (zero-time damping = off)
TDAMP = 0;
sdamp = 0.5;

% Parameters for STA/LTA Picker (TWIN=5)
BE=[0 1.0];
STA=0.007;
LTA=0.08;
TR=2;
DTR=1;
PEM=0;
PET=0;
PNL=0;
ATL=0;

% basename for first arrival picks
twin1_basename='RAJZEL_picks/Kleinneudorf_fsurf_p100/picks';

% OFFSET_MUTE
% OFFSET_MUTE = 1 mute far-offset traces
% OFFSET_MUTE = 2 mute near-offset traces
OFFSET_MUTE = 0;
OFFSETC = 50.0;

% FILTER = 0 (bandpass filter = off)
% FILTER = 1 (trapez bandpass filter = on) 
% FILTER = 2 (Butterworth low-pass filter = on)
% FILTER = 3 (Butterworth band-pass filter = on)
FILTER = 2;
order = 6;
fcmin = 40.0;
fcmax = 80.0;

% CORR = 0 (spreading correction = off)
% CORR = 1 (spreading correction (Bleistein 1986))
% CORR = 2 (spreading correction "direct wave transformation" (Forbriger et al 2014))
% CORR = 3 (spreading correction "single velocity transformation" (Forbriger et al 2014))
% CORR = 4 (hybrid transformation (Forbriger et al 2014))
% CORR = 5 (polynomial transformation)
CORR=3;
w0=1e20;        % value of w0 at w=0 to avoid singularity
fmax=120.0;
vph=140.0;
eps_corr = 1e0;
offset_corr1 = 1.0;

% load velocity picks
PICKS=0;

if(PICKS==1)
    pickdata = load('Picks_rayleigh/picks_rayleigh_1.dat');
    pick_offset = pickdata(:,1);
    pick_time = pickdata(:,2)-0.15;
end

% EST_SOURCE = 1 (estimate source wavelet in the frequency domain)
% EST_SOURCE = 2 (estimate source wavelet in the time domain)
EST_SOURCE = 0;
eps_w = 0.01;     % damping factor in Marquardt-Levenberg regularization for source wavelet estimation
eps_w_td = 100.0; 
max_offset_source = 20.0; % maximum offset of traces used for STF-inversion [m]
out_wavelet = 'wavelet/wavelet_s140';

% normalize source time function to 
% 1 - abs(maxium_value)
% 2 - sum(D_s)
NORM_SOURCE = 0;  

% damp source wavelet
DAMP_SOURCE = 0;
t0_source = 0.08;      
t1_source = 1.6;
adm_source = 1e3;  % damping parameter to smooth wavelet at the begining
adp_source = 1e3;  % damping parameter to smooth wavelet at the end

% Compare source wavelets from Matlab and DENISE
COMP_WAVELET = 0;

% Mute noisy traces during the inversion
TRACEMUTE=0;
trace_index=[1 2];

% create Near-offset section
NOS=0;

% name of input SU data file
in_data=['raw_data/SH_shots_s140/SH_s140'];

% name of input SU model file
in_model=['FD_seis/profile_s140/fwi_11_01_2019/DENISE_kleinneudorf'];

minshot = 1;          % minimum shot no. to display
maxshot = 1;          % maximum shot no. to display
maxshot_true = 50;    % true maximum number of shots

sdatname = 'source/source_kleinneudorf_s140.dat';    % specify name of source file 
rdatname = 'receiver/receiver_kleinneudorf_s140.dat';  % specify name of receiver file 

if(WIG==6)
    
    % define maximum, minimum and increment p-values
    pmax = 1.0./20.0;
    pmin = 1.0./500.0;
    np = 10000;         % number of p-values to evaluate
    mu=0.03;           % regularization parameter
    sol='adj'; 
    pfmin = 1.0;
    pfmax = 200.0;
    tau_max=0.12;

    % define p-values to compute
    np = np - 1;
    dp = (pmax-pmin)./np; % p-value increment
    p=pmin:dp:pmax;
    
end

% open FD spike wavelet
if(EST_SOURCE>0)

    in_spike=['FD_seis/profile_s140/Kleinneudorf_source_signal.su']; % spike wavelet    

    % load spike data
    [D_spike,Htr_spike,H_spike] = ReadSu(in_spike,'endian','l');
    
end

Fig = figure;
figure(Fig)

% loop over shots
for k=minshot:maxshot
    
  % if((k~=40)&&(k~=46)&&(k~=49))  % shots excluded from FWI
      
    if(WIG~=3)
      cla(gca);
    end
    
    % load model data
    infile_model = [in_model,'_y.su.shot',int2str(k),'.it1'];
    [D_mod,Htr,H_mod] = ReadSu(infile_model,'endian','l');
    
    % load model data
%     infile_model = [in_model,'_y.su.shot1.it1'];
%     [D_mod,Htr,H_mod] = ReadSu(infile_model,'endian','l');
    
    % load field data
    infile_data = [in_data,'_shot_',int2str(k),'.su'];
    [D_data,Htr,H_data] = ReadSu(infile_data,'endian','b');
%       [D_mod,Htr,H_mod] = ReadSu(infile_data,'endian','b');
    
    % reverse phase of D_data
    % D_data = -D_data;
    
%     tmp = D_mod(:,istart(k):iend(k));
%     clear D_mod;
%     D_mod = tmp;
%     
%     tmp = D_data(:,istart(k):iend(k));
%     clear D_data;
%     D_data = tmp;
   
    DT=H_mod(1).dt.*1e-6;      % get sample interval from model SU-Header
    DT_dat=H_data(1).dt.*1e-6; % get sample interval from data SU-Header

    [nt,ntr]=size(D_mod); % number of time samples and traces of model files
    t=1:nt;
    t=t.*DT;
    t=t';

    [nt_data,ntr_data]=size(D_data); % number of time samples and traces of data files
    t1=1:nt_data;
    t1=t1.*DT_dat;
    t1=t1';
    
    % calculate source coordinates
    x0 = 2.25;
    z0 = 89.74 + 6.11;
    xprofil_shift = 0.0;   
    
    xsrc(k) = Htr(1).SourceX / 100.0 + x0;
    zsrc(k) = -Htr(1).SourceSurfaceElevation / 100.0 + z0;
    
    for i=1:ntr_data
        
        xrec(i) = Htr(i).GroupX / 100.0 + x0;
        zrec(i) = -Htr(i).ReceiverGroupElevation / 100.0 + z0;           
        
    end
       
    
    offset = xsrc(k) - xrec;
      
    if(WRITE_ACQ==1)
        
        geoph1=[xrec' zrec'];
        dlmwrite(['receiver/receiver_kleinneudorf_s140_shot_',int2str(k),'.dat'],geoph1,'delimiter','\t','precision','%.3f');
        
    end
    
    %clear xrec yrec;
    
    % Interpolate data DT_dat -> DT_DENISE
    % -----------------------------------
    
    for i=1:ntr
         D_p_i(:,i) = interp1(t1,D_data(:,i),t,'spline'); 
    end
  
    clear D_data;
    D_data = D_p_i(1:nt,:);
    clear D_p_i;
    
    % apply time delay
    if(DELAY==1)
      D_data = delay_seis(delay,DT,ntr,nt,D_data,gamma_delay);
      %[D_mod] = delay_seis(nt.*DT./100,DT,ntr,nt,D_mod,gamma_delay);
    end
    
    % apply static corrections for receivers
    if(STAT_REC==1)
      D_data = stat_corr(rec_stat,DT,ntr,nt,D_data,1);
    end
    
    % apply static corrections for sources
    if(STAT_SRC==1)
      D_data = stat_corr(src_stat(k),DT,ntr,nt,D_data,2);
    end
    
    % apply bandpass filter
    % ---------------------
    if(FILTER==1)  % simple trapez filter
        D_data =  bp_filter(D_data,DT,f1,f2,f3,f4);
    end    
    
    if(FILTER==2 || FILTER==3)  % Butterworth filter
        
        % Butterworth low-pass filter
        D_data =  butter_filter(ntr,D_data,DT,fcmax,order,1);
        D_mod =  butter_filter(ntr,D_mod,DT,fcmax,order,1);
        
        if(FILTER==3)
            % Butterworth high-pass filter
            D_data =  butter_filter(ntr,D_data,DT,fcmin,order,2);
            D_mod =  butter_filter(ntr,D_mod,DT,fcmin,order,2);
        end        
               
    end
    
    % apply geometrical spreading correction to field data
    D_data = spread_corr(ntr, nt, DT, fmax, w0, offset, t, D_data, D_mod, CORR, eps_corr, vph, offset_corr1);
    
    % time window for field data
    % --------------------------
    if(TWIN>0)

        % linear time window    
        if(TWIN==1)    
    
            nnncent(1) = twinc;
            nnnorm(1) = twinc+twinp;
            nnnorm0(1) = twinc-twinm;

            for i=1:ntr
                nnnorm0(i)  =  twinc + (offset(i))./vnmo;
            end

            % calculate discrete time samples
            nncent = round(nnncent./DT);
    
        end

        % read first arrival picks
        if((TWIN==2)||(TWIN==3))
    
          % load picks of first arrival for each shot
          tmp = load([twin1_basename,'_',int2str(k),'.dat']);

          % offset_filt_m = DH .* (isrc(k) - tmp(:,1));
%          time_filt_m = DT .* tmp(:,1);
          nnnorm0 = tmp(:,1) + delay; 
          nnnorm1 = tmp(:,1) + delay + twinp; 
          % clear tmp;

          % use spline interpolation to estimate time window values 
          % at each trace
%           nnnorm0 = spline(offset_filt_m,time_filt_m,offset);
%    
%           nnnorm0 = nnnorm0 + delay_tw;
          
          if(TWIN==3)
              nnnorm1(:) = twinp_const;
          end
          
        end
        
        if(TWIN==5)
            for i=1:ntr
                sig = D_data(:,i);
                [trigt,stav,ltav,ratio,tim1,tim2,tim3,trigs,dtrigs]=stalta(sig,DT,BE,STA,LTA,TR,DTR,PEM,PET,PNL,ATL);
                nnnorm1(i) = trigs(1,1).*DT + twinp;
                nnnorm0(i) = -1.0;
            end
        end

        % calculate discrete time samples
        nnorm0 = round(nnnorm0./DT);
        nnorm1 = round(nnnorm1./DT);

    end
    
    % apply time damping
    % ------------------
    if(TDAMP>0)

        for i=1:ntr  
            D_data(:,i) = D_data(:,i).*exp(-sdamp.*t(:));
            D_mod(:,i) = D_mod(:,i).*exp(-sdamp.*t(:));
        end

    end

    % Data Normalization
    if(NORM==1) % Normalization to maximum value of each trace
        for i=1:ntr 
                 
            max_data(i) = max(abs(D_data(:,i)));
            max_FD(i) = max(abs(D_mod(:,i)));
         
            % scale data amplitudes
            D_data(:,i) = D_data(:,i) .* (max_FD(i)./max_data(i));
         
        end
    end
    
    if(NORM==2) % Normalization to maximum value of each shot
        
            max_data = max(max(abs(D_data)));
            max_FD = max(max(abs(D_mod)));
         
            % scale data amplitudes
            D_data = D_data ./ max_data;
            D_mod = D_mod ./ max_FD;
         
    end
    
    if(TRACEMUTE==1)
    
        trace_index_tmp = trace_index;
        
        % mute noisy traces
%         if(k==12)
%            trace_index_tmp = [1];
%         end
%         
%         if(k==13)
%            trace_index_tmp = [51 52 60];
%         end
%         
%         if(k==15)
%            trace_index_tmp = [65];
%         end
%         
%         if(k==16)
%            trace_index_tmp = [58 59 60];
%         end
%         
%         if(k==19)
%            trace_index_tmp = [72];
%         end
%         
%         if(k==24)
%            trace_index_tmp = [55 59 66];
%         end
%         
%         if(k==25)
%            trace_index_tmp = [68 71];
%         end
%         
%         if(k==40)
%            trace_index_tmp = [18];
%         end
%         
%         if(k==41)
%            trace_index_tmp = [12];
%         end
%         
%         if(k==43)
%            trace_index_tmp = [12 13];
%         end
%         
%         if(k==46)
%            trace_index_tmp = [10 12];
%         end
%         
%         if(k==47)
%            trace_index_tmp = [10 11 12 14 16];
%         end
        
        nkill = length(trace_index_tmp);
    
        color_trace(1:ntr) = 'k';
        
         for i=1:nkill
%             D_data(:,trace_index_tmp(i))=0.0;
%             D_mod(:,trace_index_tmp(i))=0.0;
             color_trace(trace_index_tmp(i)) = 'r';

         end


    
    end
    
    % mute far-offset 
    if(OFFSET_MUTE==1)
        for i=1:ntr
            if(abs(offset(i))>=OFFSETC)
                D_data(:,i)=0.0;
                D_mod(:,i)=0.0;
            end
        end
    end
    
    % mute far-offset 
    if(OFFSET_MUTE==2)
        for i=1:ntr
            if(abs(offset(i))<OFFSETC)
                D_data(:,i)=0.0;
                D_mod(:,i)=0.0;
            end
        end
    end
    
    % gain model and field data
    if(GAINDATA==1)
        gain=exp(again.*t);
        for i=1:ntr
            D_data(:,i)=gain.*D_data(:,i);
            D_mod(:,i)=gain.*D_mod(:,i);
        end
    end
   
    % estimate source wavelet (frequency domain)
    if(EST_SOURCE==1)
      [D_mod,D_s] = Est_source(nt,ntr,D_mod,D_data,D_spike,NORM_SOURCE,DAMP_SOURCE,out_wavelet,k,eps_w,max_offset_source,offset);
    end
   
    % estimate source wavelet (time domain)
    if(EST_SOURCE==2)
        
        for i=1:ntr
                [a(:,i),lags] = xcorr(D_data(:,i),D_mod(:,i),nt);
                [b(:,i),lags1] = xcorr(D_mod(:,i),D_mod(:,i),nt);
        end
        
        nc=length(a);
        
        for j=1:nc
            sumen(j)=0.0;
            sumed(j)=0.0;
        end
        
        % calculate nominator and denominator for Wiener filter
        Ebar = 0.0;
        for i=1:ntr    
            for j=1:nc
                
                if(abs(offset(i))<=max_offset_source)
                    sumen(j) = sumen(j) + a(j,i);
                    sumed(j) = sumed(j) + b(j,i);
                    Ebar = Ebar + abs(sumed(j));
                end
         
            end
        end
        
        Ebar = Ebar./(nc.*ntr);
        
        % calculate new source wavelet
        D_s = (ntr.*sumen)./((ntr.*sumed)+(eps_w.*ntr.*Ebar));
        D_s(1:round(length(D_s)./2))=[];
%         D_s1 = conv(D_s,D_spike,'same');
%         clear D_s;
%         D_s = D_s1;
        
        % calculate seismogram with new source wavelet
        D_s_fft = fft(D_s');
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
    
    
    % calculate seismogram with new source wavelet
%     in_spike=['Cachtice_data/profile_1/Untergrundmodell_source_signal.0.su.shot1']; % spike wavelet    
% 
%     % load spike data
%     [D_spike,H_spike] = readsegy(in_spike);
%     D_spike_fft = fft(D_spike);
%         
%     for i=1:ntr
%         D_model_fft = fft(D_mod(:,i));
%         D_sp = D_model_fft.*D_spike_fft;
%         D_mod(:,i) = real(ifft(D_sp));
%     end
    
    
    % Data Normalization
    if(NORM==1) % Normalization to maximum value of each trace
        for i=1:ntr 
                 
            max_data(i) = max(abs(D_data(:,i)));
            max_FD(i) = max(abs(D_mod(:,i)));
         
            % scale data amplitudes
            D_data(:,i) = D_data(:,i) .* max_FD(i) ./ max_data(i);
         
        end
    end
    
    if(TWIN)
        
                % apply time damping
        for i=1:ntr
     
            if(nnorm0(i)<=0)
                nnorm0(i)=1;
            end
            
             % apply time window before first arrival to data for noise
             % suppression
             if(nnorm0(i)<nt)                              
                D_data(1:nnorm0(i),i) = D_data(1:nnorm0(i),i).*exp(-adm.*(t(1:nnorm0(i))-t(nnorm0(i))).^2);
                D_mod(1:nnorm0(i),i) = D_mod(1:nnorm0(i),i).*exp(-adm.*(t(1:nnorm0(i))-t(nnorm0(i))).^2);
             end                         
             
             % apply time window after data traces
             if(nnorm1(i)>0)
                D_data(nnorm1(i):nt,i) = D_data(nnorm1(i):nt,i).*exp(-adp.*(t(nnorm1(i):nt)-t(nnorm1(i))).^2);
                D_mod(nnorm1(i):nt,i) = D_mod(nnorm1(i):nt,i).*exp(-adp.*(t(nnorm1(i):nt)-t(nnorm1(i))).^2);
             end
     
        end
        
    end
    
    if(FK_FILT==1)
        [S,kn,f] = fk_filter(D_data,DT,DX,Lham,vmin,vmax);
    end
    
    if(ENVELOPE==1)
        [D_data] = envelope_data(D_data,ntr,nt);
    end
    
    if(DAMP_SW==1)
        %[D_data] = damp_sw(D_data,offset,t,maxoffset_sw,vmin_sw,vmax_sw,gamma_sw);
        [D_data] = damp_sw(D_data,offset,t,maxoffset_sw,vmin_sw1,vmax_sw1,gamma_sw);
        [D_mod] = damp_sw(D_mod,offset,t,maxoffset_sw,vmin_sw1,vmax_sw1,gamma_sw);
    end
    
    if(NORM==2) % Normalization to maximum value of each shot
        
            max_data = max(max(abs(D_data)));
            max_mod = max(max(abs(D_mod)));
         
            % scale data amplitudes
            D_data = D_data ./ max_data;
            D_mod = D_mod ./ max_mod;
            
    end
    
    % Data Normalization
    if(NORM==3) % Normalization to maximum value of each trace of each dataset
        for i=1:ntr 
                 
            max_data(i) = max(abs(D_data(:,i)));
            max_FD(i) = max(abs(D_mod(:,i)));
         
            % scale data amplitudes
            D_data(:,i) = D_data(:,i) ./ max_data(i);
            D_mod(:,i) = D_mod(:,i) ./ max_FD(i);
         
        end
    end
    
    
    % extract and write frequency data for each shot in separate binary file
    if(WRITE_DATA==1)  

        for k1=1:length(f)
            
            imfile1=['DFT/descramble_p_shot_',int2str(k),'_stage_1_nfreq_',int2str(k1),'.bin'];
            
            h1 = 1;
            % extract FD data by DFT
            for i=1:ntr
                
                dftr = 0.0;
                dfti = 0.0;
                for j=1:nt
                    
                    t = j * DT;

                    trig1 = cos(2.0*pi*t*f(k1));
                    trig2 = sin(2.0*pi*t*f(k1));
                    
                    dftr = dftr + trig1 .* D_data(j,i);
                    dfti = dfti + trig2 .* D_data(j,i);
                    
                end
                
                DFT(h1) = dftr;
                h1 = h1 + 1;
                DFT(h1) = dfti;
                h1 = h1 + 1;
              
            end
            
            % output of DFT results
            fid1=fopen(imfile1,'w','ieee-le');
            fwrite(fid1,DFT,'float')
            fclose(fid1);
            
            clear DFT;
  
        end
        
    end
    
    % write data for each shot in separate SU file
    if(WRITE_DATA==2)  

      imfile1=[out_su_file,'_y.su.shot',int2str(k)];
      [H_shot] = make_su_file(imfile1,D_data,DT,offset.*(abs(SCALCO)),SCALCO);

    end
    
    % plot data of individual shots
    % -----------------------------

    if(MISFIT==2)
        D_res = D_mod-D_data; 
    end
    
    if(MISFIT==5)
        
        for i=1:ntr
           
           abs_D_mod = sqrt(D_mod(:,i)'*D_mod(:,i));
           abs_D_data = sqrt(D_data(:,i)'*D_data(:,i));
           dotp_mod_data = D_mod(:,i)'*D_data(:,i);
            
           for j=1:nt
              D_res(j,i) = (D_mod(j,i)./abs_D_mod.^2).*(dotp_mod_data./(abs_D_mod.*abs_D_data))-(D_data(j,i)./(abs_D_mod.*abs_D_data)); 
           end
           
        end
        
    end
    
    if(WIG==1)
       
        subplot 121
       
        if(NORM~=0)
            for i=1:dtr:ntr 

                plot(xrec(i)+(D_data((1:dresamp_nt:nt),i)./(clip.*max(abs(D_data((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),'k-','Linewidth',3.0);
                hold on;
                plot(xrec(i)+(D_mod((1:dresamp_nt:nt),i)./(clip.*max(abs(D_mod((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),'r-','Linewidth',3.0);

            end
        end
        
        if(NORM==0)
            for i=1:dtr:ntr 

                plot(xrec(i)+(D_data((1:dresamp_nt:nt),i)),t(1:dresamp_nt:nt),'k-','Linewidth',3.0);
                hold on;
                plot(xrec(i)+D_mod((1:dresamp_nt:nt),i),t(1:dresamp_nt:nt),'r-','Linewidth',3.0);

            end
        end
        
        if(TWIN>0)
            hold on;
            plot(xrec,nnnorm0,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);
            hold off;
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
           axis tight
           
           if(Tmax>0.0)
               axis([min(xrec) max(xrec) 0.0 Tmax]);
           end
           
           xlabel('offset [m]');
           ylabel('time [s]');
           iter_text=['Data comp. shot no. ',int2str(k)];
           title(iter_text);
           legend('field data','model data',1);
           
           hold off;
        
        subplot 122
        
        for i=1:dtr:ntr 
                 
            plot(offset(i)+D_res(:,i)./(clip.*max(abs(D_data(:,i)))),t,'k-','Linewidth',3.0);
            hold on;
       
        end
        
%          if(TWIN>0)
%             hold on;
%             plot(offset,nnnorm0,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);
%             hold off;
%          end
        
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
           axis tight

           if(Tmax>0.0)
               axis([min(offset) max(offset) 0.0 Tmax]);
           end
           
           xlabel('offset [m]');
           ylabel('time [s]');
           iter_text='Data residuals';
           title(iter_text);
        
           hold off;
    end

    % Display data as image plot (field data vs. model data)
    if(WIG==0)
        
         %colormap(gray);
         colormap(seismic);
         imagesc(xrec,t,D_data);
         %imagesc(offset,t,D_data);
         % imagesc(D_data);
         caxis([caxis_value1 caxis_value2]);
         
         if(TWIN>0)
            hold on;
            plot(xrec(2:ntr-1),nnnorm0(2:ntr-1),'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
            plot(xrec(2:ntr-1),nnnorm1(2:ntr-1),'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
            hold off;
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
           
           if(Tmax>0.0)
               axis([min(xrec) max(xrec) Tmin Tmax]);
           end

           xlabel('x [m]');
           ylabel('time [s]');
           iter_text=['Field data shot no. ',int2str(k)];
           title(iter_text);
         
    end
    
    % Display data as image plot (field data + model data)
    if(WIG==2)
        
         subplot 121
         %colormap(gray);
         colormap(seismic);
         imagesc(xrec,t,D_data);
         caxis([caxis_value1 caxis_value2]);
         
         if(TWIN>0)
            hold on;
            plot(xrec,nnnorm0,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
            plot(xrec,nnnorm1,'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
            hold off;
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

           if(Tmax>0.0)
               axis([min(xrec) max(xrec) Tmin Tmax]);
           end
           
           xlabel('x [m]');
           ylabel('time [s]');
           iter_text=['Field data shot no. ',int2str(k)];
           title(iter_text);
         
         subplot 122
         
         %colormap(gray);
         colormap(seismic);
         imagesc(xrec,t,D_mod);
         caxis([caxis_value1 caxis_value2]);
         if(TWIN>0)
            hold on;
            plot(xrec,nnnorm0,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
            plot(xrec,nnnorm1,'ro','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
            hold off;
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

           if(Tmax>0.0)
               axis([min(xrec) max(xrec) Tmin Tmax]);
           end
           
           xlabel('x [m]');
           ylabel('time [s]');
           iter_text=['Model data shot no. ',int2str(k)];
           title(iter_text);
         
    end
    
    % Plot AVO for each shot
    if(WIG==3)
        
        % calculate RMS amplitudes 
        for i=1:ntr
            AVO_data(i) = sqrt(sum(D_data(:,i).^2)./nt); 
            AVO_mod(i) = sqrt(sum(D_mod(:,i).^2)./nt); 
        end
        
        plot(offset(1:ntr-1),AVO_data(1:ntr-1),'r-o','Linewidth',4.0,'MarkerFaceColor','r');
        hold on;
        plot(offset(1:ntr-1),AVO_mod(1:ntr-1),'b-o','Linewidth',4.0,'MarkerFaceColor','b');
        
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
        %axis ij
        axis tight

        xlabel('Offset [m]');
        ylabel('RMS amplitude []');
        iter_text=['AVO for shot no. ',int2str(minshot),' - ',int2str(maxshot)];
        title(iter_text);
        legend('field data','model data',1);
        
        % write picks of first arrivals (- picks) and maximum time window values (+ picks) to ASCII file   
        imfile=['AVO_shot_',int2str(k),'.dat'];

        fid = fopen(imfile,'w');
            for i=1:ntr
                fprintf(fid,'%e \t %e \n',offset(i),AVO_mod(i));
            end
        fclose(fid);
        
    end
    
    % Display data as image plot (field data only)
    if(WIG==4)
        
         %colormap(gray);
         %colormap(seismic);
         %imagesc(offset,t,D_data);
         
%          if(NORM==2)
%              clip=0;
%          end
         
         for i=1:dtr:ntr 
            hold on;
            
            if(TRACEMUTE==0)
               plot(xrec(i)+(D_data((1:dresamp_nt:nt),i)./(clip.*max(abs(D_data((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),'k-','Linewidth',3.0);                
            end
             
            if(TRACEMUTE==1)
               plot(xrec(i)+(D_data((1:dresamp_nt:nt),i)./(clip.*max(abs(D_data((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),color_trace(i),'Linewidth',3.0);
               %plot(offset(i)+(D_data((1:dresamp_nt:nt),i)./(clip.*max(abs(D_data((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),color_trace(i),'Linewidth',3.0);
            end

         end                 
         
         if(TWIN>0)
            hold on;
            plot(offset,nnnorm0,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);
            plot(offset,nnnorm1,'bo','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
            hold off;
         end
         
         if(PICKS==1)
             plot(pick_offset,pick_time,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);
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
           
           if(Tmax>0.0)
             axis([min(xrec) max(xrec) 0.0 Tmax]);
             %axis([min(offset) max(offset) 0.0 Tmax]);
           end

           xlabel('xrec [m]');
           ylabel('time [s]');
           iter_text=['Field data shot no. ',int2str(k)];
           title(iter_text);
         
    end
    
    % Calculate and plot phase-velocity frequency spectra
    if(WIG==6)
        
       dx = offset(2) - offset(1);
       [S_data,f] = Vf_spectra(D_data,DT,dx,200.0,offset,cphase);
       [S_mod,f] = Vf_spectra(D_mod,DT,dx,200.0,offset,cphase);
       
       colormap(jet(256));
       
       subplot(1,2,1)
       imagesc(f,cphase,abs(S_data'));
       %axis([pmin pmax 0 tau_max]);
       hold on;
       
       set(gca,'YDir','normal');
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
       %axis ij

       ylabel('Phase velocity [m/s]');
       xlabel('Frequency [Hz]');
       iter_text=['Field data shot no. ',int2str(k)];
       title(iter_text);
       
       subplot(1,2,2)
       imagesc(f,cphase,abs(S_mod'));
       %axis([pmin pmax 0 tau_max]);
       hold on;
       
       set(gca,'YDir','normal');
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
       %axis ij

       ylabel('Phase velocity [m/s]');
       xlabel('Frequency [Hz]');
       iter_text=['Modelled data shot no. ',int2str(k)];
       title(iter_text);
         
    end
   
% Estimate Qs-values by spectral division
    if(WIG==7)
 
       [f,Q,U,As] = Q_spec_div(D_data,D_s,DT,100.0,offset);
       
       colormap(jet(256));
      
       [n1,n2]=size(Q);
       clear n1;
       
       for i=1:n2
          plot(f,Q(:,i),'r-','Linewidth',2.0,'MarkerFaceColor','r');
          hold on;
       end
       
       set(gca,'YDir','normal');
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
       %axis ij

       ylabel('Qs []');
       xlabel('Frequency [Hz]');
       iter_text=['Field data shot no. ',int2str(k)];
       title(iter_text);
                
    end
    
    % Estimate data fit in the frequency domain
    if(WIG==8)
       
        for i=1:ntr
            
           nf = 4*2^nextpow2(length(D_data(:,i)));
           n2 = nf/2+1;
           X = fft(D_data(:,i),nf); 
           X1 = fft(D_mod(:,i),nf); 
           fa = [0:1:nf-1]'/nf/DT;

           Phase_Shift = exp(-i*2*pi*fa*DT);
           X = X.*Phase_Shift;
           X1 = X1.*Phase_Shift;
           n2 = floor( fmax* (nf*DT)) +1;
           X = X(1:n2,1);
           X1 = X1(1:n2,1);
       
           Theta = angle(X.*conj(X1));
           Theta(1,1) = 0.;
           Theta = Theta*180/pi;
           phase = Theta;
           
           amp_res(:,i) = abs(X-X1);
           phase_res(:,i) = phase;
           
           clear X fa Theta;
        end
       
       f = [1:1:n2]'; f = (f-1)/nf/DT;
        
       colormap(jet(256));
       
       subplot(1,2,1)
       imagesc(offset,f,amp_res);
       cbar_handle=colorbar('SouthOutside');
       set(get(cbar_handle,'xlabel'),'string','Amplitude-error []','fontsize',FSize,'FontWeight','bold');
       hold on;
       
       set(gca,'YDir','normal');
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
       %axis ij

       xlabel('Offset [m]');
       ylabel('Frequency [Hz]');
       iter_text=['Amplitude residuals shot no. ',int2str(k)];
       title(iter_text);
       
       subplot(1,2,2)
       imagesc(offset,f,phase_res);
       caxis([-180,180]);
       cbar_handle=colorbar('SouthOutside');
       set(get(cbar_handle,'xlabel'),'string','Phase-error [�]','fontsize',FSize,'FontWeight','bold')
       hold on;
       
       set(gca,'YDir','normal');
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
       %axis ij

       xlabel('Offset [m]');
       ylabel('Frequency [Hz]');
       iter_text=['Phase residuals shot no. ',int2str(k)];
       title(iter_text);
       
    end
    
    % Display fk-spectrum (field data only)
    if(WIG==9)
        
         %colormap(gray);
         colormap(seismic);
         if(FK_FILT==0)
            [S,kn,f] = fk_spectra(D_data,DT,DX,Lham);
         end
         
         imagesc(offset,t,real(S(1:nt,1:ntr)));
         caxis([caxis_value1 caxis_value2]);
         
        
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
           %axis([min(kn) max(kn) 0 fmaxfk]);
           axis ij

           if(Tmax>0.0)
               axis([min(offset) max(offset) 0.0 Tmax]);
           end

           xlabel('offset [m]');
           ylabel('time [s]');
           iter_text=['Field data shot no. ',int2str(k)];
           title(iter_text);
           
           clear kn f;
         
    end
    
    % Plot acquisition geometry
    if(WIG==10)       
         
         plot(xsrc(k),ysrc(k),'ro');
         hold on;
         plot(xrec,yrec,'ks');
         plot(0,0,'gs');
        
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
           %axis([min(kn) max(kn) 0 fmaxfk]);
           axis ij

           xlabel('x [m]');
           ylabel('y [m]');
           iter_text=['Acquisition geometry shot no. ',int2str(k)];
           title(iter_text);
         
    end
    
    % Display stacked data
    if(WIG==11)
        
         if(k==1)
             stack = D_data;
         end
        
         if(k~=14)
             stack = stack + D_data;
         end
        
         colormap(gray);
         %colormap(seismic);
         imagesc(xrec,t,stack);
         caxis([caxis_value1 caxis_value2]);
              
         
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
           
           if(Tmax>0.0)
             axis([min(xrec) max(xrec) 0.0 Tmax]);
           end

           xlabel('offset [m]');
           ylabel('time [s]');
           iter_text=['Field data shot no. ',int2str(k)];
           title(iter_text);
         
    end
    
    if(WIG==12)
       
        if(NORM~=0)
            for i=1:dtr:ntr 

                plot(xrec(i)+(D_data((1:dresamp_nt:nt),i)./(clip.*max(abs(D_data((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),'k-','Linewidth',3.0);
                hold on;
                plot(xrec(i)+(D_mod((1:dresamp_nt:nt),i)./(clip.*max(abs(D_mod((1:dresamp_nt:nt),i))))),t(1:dresamp_nt:nt),'r-','Linewidth',3.0);

            end
        end
        
        if(NORM==0)
            for i=1:dtr:ntr 

                plot(xrec(i)+(D_data((1:dresamp_nt:nt),i)),t(1:dresamp_nt:nt),'k-','Linewidth',3.0);
                hold on;
                plot(xrec(i)+D_mod((1:dresamp_nt:nt),i),t(1:dresamp_nt:nt),'r-','Linewidth',3.0);

            end
        end
        
        if(TWIN>0)
            hold on;
            plot(xrec,nnnorm0,'go','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);
            hold off;
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
           axis tight
           
           if(Tmax>0.0)
               axis([min(xrec) max(xrec) 0.0 Tmax]);
           end
           
           xlabel('offset [m]');
           ylabel('time [s]');
           iter_text=['Data comp. shot no. ',int2str(k)];
           title(iter_text);
           legend('field data','model data',1); 
           hold off;
    end
    
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperSize', [5 7]);

      set(Fig,'position',[0 0, screen_x screen_z])
      set(Fig,'PaperPositionMode','Auto')       

      output=['pics/Kleinneudorf_data_shot_',int2str(k)];
      %saveas(Fig,output,'psc2'); 
      saveas(Fig,output,'png');
    
   if(TWIN>0)
        
        % write picks of first arrivals (- picks) and maximum time window values (+ picks) to ASCII file   
        imfile=['picks/picks',int2str(k),'.dat'];

        fid = fopen(imfile,'w');
            for i=1:ntr
                fprintf(fid,'%e \n',nnnorm1(i));
            end
        fclose(fid);
    
   end 
   
    % create near offset section
    % --------------------------
    if(NOS==1)
        
        if(k<=(maxshot_true-2))
          section_NOS(:,k) = D_data(:,k+2);    
        end
        
        if(k>=(maxshot_true-2))
          section_NOS(:,k) = D_data(:,k-2);    
        end
        
    end
    
   % calculate phase error 
   % ---------------------
   if(FWI_phase_error==1)
       for i=1:ntr
         phase_diff(i,k)  = calc_phase_error(DT,D_mod(:,i),D_data(:,i),DT,80.0,1,fphase);
       end
   end
    
    clear offset;
    
    clear D_data;
    clear D_mod;
    
    if(WIG~=3)
        hold off;
    end
    
    if(TRACEMUTE==1)
    
        for i=1:ntr
        
            trace_mute(i)=0;
        
                for j=1:nkill
                    if(i==trace_index_tmp(j))
                        trace_mute(i)=1;
                    end
                end
        
        end
    
        % write tracekill data to ASCII file   
        trace_mute_file=['trace_kill.dat'];

        if(k==minshot)
            fid = fopen(trace_mute_file,'w');
        else
            fid = fopen(trace_mute_file,'a');
        end
    
        for i=1:ntr
             fprintf(fid,'%d \t',trace_mute(i));
        end
        
        fprintf(fid,'\n');
      
        fclose(fid);
    
        clear trace_index_tmp;
        
    end
    
    pause(0.2);
    
    clear xrec zrec noffset offset geoph1;
    
    fclose('all');
    
  % end % end of if-statement excluding neglected shots
    
end

if(WRITE_ACQ==1)
    
    imfile=['source/source_kleinneudorf_s140.dat'];
    time_shift = 9e-3;
    ntr = length(xsrc);
    fid = fopen(imfile,'w');
    fprintf(fid,'%d \n',ntr);
    for i=1:ntr
    fprintf(fid,'%6.2f \t %6.2f \t %6.2f \t %6.3f \t %6.2f \t %6.2f \t %6.2f \t %d \n',xsrc(i),0.0,zsrc(i),time_shift,10.0,1.0,0.0,3);
    end

    fclose(fid);
    
end

if(COMP_WAVELET==1)
    
    Fig = figure;
    figure(Fig)
    
    data = load('Ainos_data/wavelet_ainos_shot_1.dat');
    plot(t,D_s./max(abs(D_s)),t,data,'Linewidth',2);

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
    %axis ij
    axis tight

    xlabel('Time [s]');
    ylabel('Amplitude []');
    iter_text=['Estimated source wavelets'];
    title(iter_text);
    legend('Matlab','DENISE',1);

    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [5 7]);

    set(Fig,'position',[0 0, screen_x screen_z])
    set(Fig,'PaperPositionMode','Auto')       

    output=['pics/wavelet_Matlab_DENISE_comp'];
    %saveas(Fig,output,'psc2'); 
    %saveas(Fig,output,'png');
    
end

if(NOS==1)

   Fig = figure;
   figure(Fig)

   % plot near offset section
   % ------------------------
   imagesc(section_NOS);
   caxis([caxis_value1 caxis_value2]);
   
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

   %xlabel('offset [m]');
   %ylabel('time [s]');
   iter_text=['Near-Offset Section '];
   title(iter_text);
            
   set(gcf, 'PaperUnits', 'inches');
   set(gcf, 'PaperSize', [5 7]);

   set(Fig,'position',[0 0, screen_x screen_z])
   set(Fig,'PaperPositionMode','Auto')       

   output=['pics/NOS'];
   %saveas(Fig,output,'psc2'); 
   %saveas(Fig,output,'png');       

end

% plot phase error
if(FWI_phase_error==1)

   % check for NaN and set phase error to 0.0
   for i=1:(maxshot-minshot)
        for j=1:ntr
            if(isnan(phase_diff(j,i))==1)
                phase_diff(j,i) = 0.0;
            end
        end
   end 
    
   Fig = figure;
   figure(Fig)

   % plot near offset section
   % ------------------------
   
   imagesc(phase_diff);
   caxis([-180.0 180.0]);
   colormap(seismic);
   cbar_handle=colorbar('EastOutside');
   set(get(cbar_handle,'ylabel'),'string','Phase error [�]','fontsize',FSize,'FontWeight','bold');
   
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
   axis ij;
   axis equal;
   axis tight;

   xlabel('source #');
   ylabel('receiver #');
   iter_text=['Phase error'];
   title(iter_text);
            
   set(gcf, 'PaperUnits', 'inches');
   set(gcf, 'PaperSize', [5 7]);

   set(Fig,'position',[0 0, screen_x screen_z])
   set(Fig,'PaperPositionMode','Auto')       

   output=['pics/phase_error'];
   % saveas(Fig,output,'psc2'); 
   %saveas(Fig,output,'png');       

end



