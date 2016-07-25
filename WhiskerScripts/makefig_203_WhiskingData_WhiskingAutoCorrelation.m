% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_201 and generates plots that are
% shown in Figure 6G of the manuscript: autocorrelation of whisking stimulus
% during whisking (high amplitude) and at all times.


clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;
Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ; % Cell indices (used in file names) 
lgnd = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

sr = 1000 ; % (Hz) sampling rate

for i = 1:7
    load([datadir 'VPM_cell_' Names{i} '_autocorrelation.mat']) ;

    figure ;
    plot(tau,mean(ACall,2),'k',tau,mean(ACisw,2),'r') ;
    legend('all','whisking') ;
    xlabel('time lag [s]') ;
    ylabel('autocorrelation') ;
    title(['cell ' lgnd{i}]) ;
end

