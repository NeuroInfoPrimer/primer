% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_203 and generates plots that are
% shown in Figure 7A of the manuscript: The spike triggered average
% feature.

% In addition plotted are the stimulus distribution projected on the STA,
% the spike conditioned stimulus distribution and the probability of a
% spike (per time bin) given a stimulus

clear ;
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;
sr = 1000 ;
N = 150 ;
ds = 2 ;
dt = ds/sr ;

fntsz = 10 ;
iJK = 5 ;
for i = 1:7
    load([datadir 'VPM_cell_' Names{i} '_sta_JK_' num2str(iJK) '.mat']) ; 
    
    subplot(2,7,i) ;
    
    plot(ctr_plot,Pf_plot,'k','LineWidth',2) ; hold on ; 
    plot(ctr_plot,Pfs_plot,'r','LineWidth',2) ; hold on ;
    plot(ctr_plot,sta_model_rect_norm(ctr_plot)*max(Pf)/max(Psf),'b','LineWidth',2) ;
    if i == 1
        legend('P(z)','P(z|spike)','P(spike|z)') ;
    end
    plot(ctr,Psf/max(Psf)*max(Pf),'ob','MarkerSize',3) ; hold on ;
    set(gca,'FontSize',fntsz) ; axis tight square ; box off 
    ylim([0 1.2]) ;
    title(lgnd{i},'FontSize',fntsz) ;

    subplot(2,7,i+7) ;
    plot(tp,sta,'k','LineWidth',3) ; 
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ; axis tight square ; box off 
    if i == 1
        legend('STA') ;
    end
    xlabel('t (ms)','FontSize',fntsz) ;
    ylim([-0.2 0.2]) ;
    xlim([-N*dt 0]) ;
end