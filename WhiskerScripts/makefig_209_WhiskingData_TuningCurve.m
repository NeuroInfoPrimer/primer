% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_209 and generates the plots
% showing the whisker phase tuning curves of all 7 thalamic cells. examples
% are shown in Figures 6 an 12 in the manuscript.

clear ; 

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

sr = 1000 ;               % (Hz) sampling rate
ds = 2 ;                  % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;              % delta t of stimulus

fntsz = 15 ;
p = linspace(-pi,pi,1000) ;

iJK = 5 ;
for i = 1:7
    load([datadir 'VPM_cell_' Names{i} '_tune_JK_' num2str(iJK) '.mat']) ;
    subplot(1,7,i) ;
    plot(p,tun_model_rect_norm(p)/dt,'k','LineWidth',3) ;
    xlim([-pi pi]) ;
    ylim([0 0.2]/dt) ;
    title(lgnd{i},'FontSize',fntsz) ;
    set(gca,'FontSize',fntsz,'xtick',-pi:(pi/2):pi,'xticklabel',{'-\pi','\pi/2','0','\pi/2','\pi'}) ;
    if i==1
        ylabel('firing rate (Hz)','FontSize',fntsz) ;
    end
end