% Build Kleinneudorf initial model with dispersion correction
% 
% Daniel Koehn
% Kiel, 23.11.2017

clear all;
close all;

FSize=25;
WRITEMODE=1;
FREE_SURF=1;
TAPER=1;

% define FD grid parameters
nx = 1952;
nz = 1392;
dh = 0.05;

sdatname = '../source/source_kleinneudorf_s140.dat';        % specify name of source file 
rdatname = '../receiver/receiver_kleinneudorf_s140.dat';    % specify name of receiver file 

% read source coordinates from source.dat
shot_data = importdata(sdatname,'\t', 1);
xshot = shot_data.data(:,1);
zshot = shot_data.data(:,3);

% read receiver positions
rec_data = load(rdatname); 
xrec = rec_data(:,1);
zrec = rec_data(:,2);

% Add additional virtual receiver to include escarpement in the model
xrec = [xrec(1:16); xrec(17)-10*dh; xrec(17:end)];
zrec = [zrec(1:16); zrec(16); zrec(17:end)];

% clean up memory
clear shot_data rec_data;

x = dh.*[1:nx];
z = dh.*[1:nz];

vs = ones(nz,nx);
rho = ones(nz,nx);

% load Vs FWI-model
nx1 = 3008;
nz1 = 1056;

file='../start/Kleinneudorf_fsurf_smooth_p500_large.vs';
disp([' loading file ' file]);
fid=fopen(file,'r','ieee-le');
vs_p500=fread(fid,[nz1,nx1],'float');
fclose(fid);

vs_p500_profile = vs_p500(65:nz1,1);
npro = length(vs_p500_profile);
vs_p500_profile(npro+1:nz) = vs_p500_profile(npro);

% define free surface via source and receiver positions
if(FREE_SURF==1)
    
    shift_fs = 4*dh;
    xfs = [xrec ; xshot];
    zfs = [zrec ; zshot];
    
    topo_geom = [xfs zfs];    
    topo_geom1 = unique(topo_geom,'rows','first');
    xyb = sortrows(topo_geom1);
    
    clear xfs zfs
    
    xfs  = xyb(:,1);
    zfs  = xyb(:,2);
    
    zfree = interp1(xfs,zfs,x,'PCHIP');
    
    plot(xfs,zfs,'bo-');
    hold on;
    plot(x,zfree,'r-');        
    axis ij;

    % define material parameters above free surface
    for i=1:nx    
        h=1;
        for j=1:nz     
            
            if(z(j)<=zfree(i)-shift_fs)
                vs(j,i) = 0.0;
                rho(j,i) = 0.0;
            end
            
            if(z(j)>zfree(i)-shift_fs)
                vs(j,i) = vs_p500_profile(h);
                h = h + 1;
            end
            
        end
    end

end

% apply dispersion correction
vs_corr = 0.0;
vs_old = vs;
vs = vs + vs_corr .* vs; 

% define density model by empirical Vs-density relation
% -----------------------------------------------------
rho_new = 0.1055 * log(vs) + 1.3871;  % density in g/cm^3
rho_new = 1000.0 * rho_new;           % density in kg/m^3

rho = rho_new;

% define taper function
if(TAPER==1)
    
    taper = vs./vs;
    
    % define material parameters above free surface
    for j=1:nz
        for i=1:nx                    
            if(vs(j,i)==0)    
                taper(j,i) = 0.0;
                rho(j,i) = 0.0;
            end
        end
    end

end

% resample vs and rho model
vs1 = vs;
% vs1 = vs(1:2:nz,1:2:nx);
% rho1 = rho(1:2:nz,1:2:nx);
% taper1 = taper(1:2:nz,1:2:nx);
% 
% nx = round(nx/2)
% nz = round(nz/2)
% dh = dh * 2
% 
% x = dh.*[1:nx];
% z = dh.*[1:nz];

figure;
colormap(jet(256));
imagesc(x,z,vs1);

hold on;

plot(xshot,zshot,'r*','Markersize',4);
plot(xrec,zrec,'w.','Markersize',6);
plot(xrec(16),zrec(16),'g.','Markersize',6);
plot(xrec(17),zrec(17),'b.','Markersize',6);

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

xlabel('Distance [m]');
ylabel('Depth [m]');
iter_text=['1D initial model'];
title(iter_text);

axis equal;
axis tight;

% write model to file
if(WRITEMODE==1)
     
  file1=['Kleinneudorf_init_s140.vs'];
  fid1=fopen(file1,'w','ieee-le');
  fwrite(fid1,vs,'float')
  fclose(fid1);
  
  file1=['Kleinneudorf_init_s140.rho'];
  fid1=fopen(file1,'w','ieee-le');
  fwrite(fid1,rho,'float')
  fclose(fid1);
  
  if(TAPER==1)
    file1=['taper_s140.bin'];
    fid1=fopen(file1,'w','ieee-le');
    fwrite(fid1,taper,'float')
    fclose(fid1);
  end
  
end