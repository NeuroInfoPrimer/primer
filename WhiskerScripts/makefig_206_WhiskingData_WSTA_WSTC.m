% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_203, 204, 205, 206 and generates
% plots that are shown in Figure 7A,C of the manuscript. 

% In the text (Fig. 7A) we show the STA and STC features computed with the
% raw stimulus (no whitening) and the features resulting from performing
% the STA and STC analysis on the whitened stimulus. To allow comparison,
% in the latter the stimulus correlations are re-introduced by multiplying
% the whitened STA and STC features by the stimulus covariance. 

% In this script, to illustrate the difference between the STA and STC 
% analysis done on the raw stimulus and the analysis done on the whitened 
% stimulus we plot the whitened STA and STC features without re-introducing
% the stimulus correlations. These features "act on" the whitened stimulus
% (i.e. the stimulus multiplied by the pseudo-square-root-inverse -- Eq. 30
% -- is projected onto them). 

% since the raw stimulus has much higher power at low frequencies, as shown
% in the principal component analysis Fig. 6A, the features computed from
% it are smooth. Whitening equalizes the power in all frequencies
% (only approximately for a finite data-set), so the STA and STC features
% plotted by this script have power also at high frequencies. 

% additionally this script plots the nonlinearity of the whitened STA and
% STC model, shown in Figure 7B of the manuscript. 

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
clr = {'r','g','m','c','k'} ;

iJK = 5 ;

for i = 1:7 
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stc_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_wstc_logl_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_wstc_L' num2str(Lopt) '_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stcsig_L' num2str(Lopt) '_JK_' num2str(iJK) '.mat']) ;
   
    zr = 3/max(bns{1}) ;
   
    tp = fliplr(tp) ;
    C = cov(S) ;
    zL = zeros(T,2) ;
    SL = S*Cpinv/N ;
    zL(:,1) = SL*staL ;
    zL(:,2) = SL*stcL(:,1) ;
    
    zx = max(abs(zL(:))) ;
    
    nbnsplot = 100 ;
    Vzplot = linspace(-zx,zx,nbnsplot) ;
    [zplot1,zplot2] = meshgrid(Vzplot,Vzplot) ;
    zplot1 = reshape(zplot1,nbnsplot^2,1) ;
    zplot2 = reshape(zplot2,nbnsplot^2,1) ;
    Mzplot = [zplot1 zplot2] ;
    stcL_surf = reshape(stc_model_rect_normL(Mzplot),nbnsplot,nbnsplot) ;
    
    figure ;
    
    imagesc(Vzplot*zr,Vzplot*zr,(stcL_surf/dt),[0 max(stcL_surf(:))/dt*1.1]) ;
    xlabel('S*staL','FontSize',fntsz) ;
    ylabel('S*stcL_1','FontSize',fntsz) ;
    axis square xy ; box on ;
    colormap('cool') ;
    xlim([-zx zx]*zr) ;
    ylim([-zx zx]*zr) ;
    colorbar('eastoutside') ;
    set(gca,'FontSize',fntsz) ;
    
    figure ;
    subplot(nsigL+1,1,1) ;
    plot(tp,staL/norm(staL),'b','LineWidth',3) ; hold on
    plot(tp,sta,'Color',[0.7 0.7 0.7],'LineWidth',3) ; hold on
    legend('staL','sta') ;
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ; axis tight square ; box off
    xlabel('t (ms)','FontSize',fntsz) ;
    ylim([-0.3 0.3]) ;
    xlim([-N*dt 0]) ;

    for iisig = 1:nsigL
        subplot(nsigL+1,1,iisig+1) ;
        plot(tp,stcL(:,iisig)/norm(stcL(:,iisig)),clr{iisig},'LineWidth',3) ; hold on
        plot(tp,vdCL(:,isigL(iisig))/norm(vdCL(:,isigL(iisig))),'--','Color',clr{iisig},'LineWidth',3) ; hold on
        if iisig <= size(stc,2)
            plot(tp,stc(:,iisig),'Color',[0.7 0.7 0.7],'LineWidth',3) ; hold on
            legend(['stcL' num2str(iisig) ' orth to staL'],['stcL' num2str(iisig)],['stc' num2str(iisig)]) ;
        else
            legend(['stcL' num2str(iisig) ' orth to staL'],['stcL' num2str(iisig)]) ;
        end
        set(gca,'FontSize',fntsz,'Yaxislocation','right') ; axis tight square ; box off
        xlabel('t (ms)','FontSize',fntsz) ;
        ylim([-0.3 0.3]) ;
        xlim([-N*dt 0]) ;
    end
    
    figure ;
    
    ctr_plot = bns_plot{1}(1:nbns_plot-1)+dbns_plot/2 ;
    ctr = bns{1}(1:nbns-1)+dbns/2 ;
    
    subplot(2,2,1) ;
    
    plot(ctr*zr,PsfL_a/dt,'b','LineWidth',2) ;
    axis tight square ; box off
    ylim([0 0.2/dt]) ;
    xlabel('S*stak','FontSize',fntsz) ;
    ylabel('firing rate [Hz]','FontSize',fntsz) ;
    
    set(gca,'FontSize',fntsz) ;
    
    subplot(2,2,3) ;
    
    plot(ctr*zr,PsfL_c/dt,'b','LineWidth',2) ;
    axis tight square ; box off
    ylim([0 0.2/dt]) ;
    xlabel('S*stck_1','FontSize',fntsz) ;
    ylabel('firing rate [Hz]','FontSize',fntsz) ;
    
    set(gca,'FontSize',fntsz) ;
    
    subplot(2,2,2) ;
    
    plot(ctr_plot*zr,PfL_a_plot/zr,'k','LineWidth',2) ; hold on ;
    plot(ctr_plot*zr,PfsL_a_plot/zr,'r','LineWidth',2) ; hold on ;
    axis tight square ; box off
    ylim([0 1.2]) ;
    xlabel('S*stak','FontSize',fntsz) ;
    ylabel('P(S*stak)','FontSize',fntsz) ;
    
    set(gca,'FontSize',fntsz) ;
    
    subplot(2,2,4) ;
    
    plot(ctr_plot*zr,PfL_c_plot/zr,'k','LineWidth',2) ; hold on ;
    plot(ctr_plot*zr,PfsL_c_plot/zr,'r','LineWidth',2) ; hold on ;
    axis tight square ; box off
    ylim([0 1.2]) ;
    xlabel('S*stck_1','FontSize',fntsz) ;
    ylabel('P(S*stck_1)','FontSize',fntsz) ;
    
    set(gca,'FontSize',fntsz) ;
    
    
end
