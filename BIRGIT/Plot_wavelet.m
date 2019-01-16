% Plot source wavelet 
%
% Daniel Koehn
% Raisdorf, 19th of July 2017

clear all;
close all;

% Locate directories for SU matlab
Dir='../SU_matlab';
path(path, strcat(Dir,'/segy'));
addpath(Dir);

% font size
FSize = 25;

% define screen size
screen_x = 1400;
screen_z = 1200;

% load wavelet
wavelet = load('wavelet/wavelet_skagerrak_shot_91.dat');

% time sample interval
DT = 0.001;
nt = length(wavelet);
t = 1:nt;
t = t .* DT;

Fig = figure;
figure(Fig)

plot(t,wavelet,'b','Linewidth',3);

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
axis tight

xlabel('Time [s]');
ylabel('Amplitude []');
iter_text='Estimated source wavelet shot no. 91';
title(iter_text);

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 7]);

set(Fig,'position',[0 0, screen_x screen_z])
set(Fig,'PaperPositionMode','Auto')       

output=['pics/Skagerrak_wavelet_shot_91'];
saveas(Fig,output,'psc2'); 

imfile1=['wavelet/wavelet_skagerrak.su'];
[H_shot] = make_su_file(imfile1,wavelet,DT,1,1);


