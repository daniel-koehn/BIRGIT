close all 
clear all 

% Q-approximation using an improved tau-method (Blanch et al., 1995) 
% otimization routine: Marquardt-Levenberg 

global L w Qf1 Qf2 

%-------------------------INPUT-PARAMETERS---------------------------- 
Q0=30.0;                   % constant Q to be approximated 
fp1=5.0;fp2=80.0;df=1.0; % within frequency range fp1,..., fp2 
% L=10;                    % number of relaxation mechanisms       
% fl_st=[0.01 0.1 1 5 10 70 200 500 1000 10000];      % L starting values for the relaxation frequencies
L=4;                    % number of relaxation mechanisms       
fl_st=[5.0 25.0 50.0 80.0]; % L starting values for the relaxation frequencies
t=2/Q0;                 % starting value for tau 
%-------------------------END: INPUT-PARAMETER----------------------- 
f=fp1:df:fp2; 
w=2*pi*f; 
Qf1=Q0+f*0;           % here constant Q-approximation, arbitrary 
                      % frequency dependency of Q possible, define 
                      % Q as function of frequency here 

%Qf1=Q0+sin(2*pi*f/fp2)*Q0; 

% defining option for otimization (see 'help foptions') 
options(1)=1; 
options(2)=0.1; 
options(3)=0.1; 
options(5)=0;     % =0 : Levenberg-Marquardt Method, =1: Gauss-Newton Method 
options(14)=500; 

x=[fl_st t]; 

% if optimization toolbox is installed use 'leastq' 
[x,options]=lsqnonlin('qstd',x,options); 

% output of results: 

Optimized_relaxation_frequencies_FL=x(1:L) 
TAU=x(L+1)
sigma=sqrt(sum((Qf2-Qf1).*(Qf2-Qf1))/size(Qf2,2))*100/Q0 
GEW_REL_STANDARDABW=sigma

% plot Q as function of frequency: 
plot(f,Qf1,f,Qf2,'r--'); 
xlabel('frequency [Hz]'); 
ylabel('quality factor (Q)');
legend('const. Qs = 10','Qs = 10 approximation');

axis([fp1 fp2 1 80]); 
