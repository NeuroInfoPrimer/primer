% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_208 and generates
% plots that are shown in Figure 9 of the manuscript: 
% panel A - the stimulus filter is compared to the STA
% panel B - the spike history filter and the exponentiated spike history
%           filter

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd  = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;
fntsz = 15 ;

iJK = 5 ;
for i = 1:7
    load([datadir 'VPM_cell_' Names{i} '_sta_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_glm_JK_' num2str(iJK) '.mat']) ;
     
    figure ;
    subplot(4,1,1) ;
    plot(flipud(tp),zeros(size(tp)),'k') ; hold on ;
    h1 = plot(flipud(tp),sta,'r','LineWidth',3) ; axis tight ; box off 
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ;
    legend(h1, 'STA') ;
    title(lgnd{i},'FontSize',fntsz) ;
    
    subplot(4,1,2) ;
    plot(tp,zeros(size(tp)),'k') ; hold on ;
    h1 = plot(tp,gg.k,'g','LineWidth',3) ; axis tight ; box off
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ;
    legend(h1, 'GLM stimulus filter') ;
    
    subplot(4,1,3) ;
    plot((tp(2)-tp(1))*gg.iht,zeros(size(gg.iht)),'k') ; hold on ;
    h1 = plot((tp(2)-tp(1))*gg.iht,(gg.ihbas*gg.ih),'b','LineWidth',3) ; axis tight ; box off
    legend(h1, 'GLM spike history filter') ;
    xlabel('time relative to spike (s)') ;
    set(gca,'FontSize',fntsz) ;

    subplot(4,1,4) ;
    h1 = plot((tp(2)-tp(1))*gg.iht,exp(gg.ihbas*gg.ih),'b','LineWidth',3) ; axis tight ; box off
    legend(h1, 'exp(GLM spike history filter)') ;
    xlabel('time relative to spike (s)') ;
    set(gca,'FontSize',fntsz) ;
end