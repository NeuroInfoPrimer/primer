% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

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
i = 3 ;
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stc_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_wstc_logl_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_wstc_L' num2str(Lopt) '_JK_' num2str(iJK) '.mat']) ;
    tp = fliplr(tp) ;
    C = cov(S) ;
    
    plot(tp,C*stcL(:,1)/norm(C*stcL(:,1)),'--r','LineWidth',3) ; hold on
    plot(tp,C*vdCL(:,isigL(1))/norm(C*vdCL(:,isigL(1))),'r','LineWidth',3) ; hold on
    
    plot(tp,C*stcL(:,2)/norm(C*stcL(:,2)),'--m','LineWidth',3) ; hold on
    plot(tp,C*vdCL(:,isigL(2))/norm(C*vdCL(:,isigL(2))),'m','LineWidth',3) ; hold on
    
    plot(tp,stc(:,1)/norm(stc(:,1)),'--k','LineWidth',3) ; hold on
    plot(tp,vdC(:,isig(1)),'k','LineWidth',3) ; hold on
    
    plot(tp,stc(:,2)/norm(stc(:,2)),'--g','LineWidth',3) ; hold on
    plot(tp,vdC(:,isig(2)),'g','LineWidth',3) ; hold on
    
%    plot(tp,stc(:,3)/norm(stc(:,3)),'--b','LineWidth',3) ; hold on
%    plot(tp,vdC(:,isig(3)),'b','LineWidth',3) ; hold on
    
    legend('whitened STC_1, prep to whitened STA',...
           'whitened STC_1',...
           'whitened STC_2, prep to whitened STA',...
           'whitened STC_2',...
           'STC_1, prep to STA',...
           'STC_1',...
           'STC_2, prep to STA',...
           'STC_2') ; 
       pL1 = (C*vdCL(:,isigL(1))/norm(C*vdCL(:,isigL(1))))'*(C*staL)/norm(C*staL) ;
       pL2 = (C*vdCL(:,isigL(2))/norm(C*vdCL(:,isigL(2))))'*(C*staL)/norm(C*staL) ;

       p1 = vdC(:,isig(1))'/norm(vdC(:,isig(1)))*(sta)/norm(sta) ;
       p2 = vdC(:,isig(2))'/norm(vdC(:,isig(2)))*(sta)/norm(sta) ;

%%
    zr = 3/max(bns{1}) ;
   
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
    
    imagesc(Vzplot*zr,Vzplot*zr,(stcL_surf/dt),[0 200]) ;
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
    plot(tp,C*staL/norm(C*staL),'b','LineWidth',3) ; hold on
    plot(tp,sta,'Color',[0.7 0.7 0.7],'LineWidth',3) ; hold on
    legend('staL','sta') ;
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ; axis tight square ; box off
    title(['pseudoinverse order: ' num2str(Lopt)],'FontSize',fntsz) ;
    xlabel('t (ms)','FontSize',fntsz) ;
    ylim([-0.3 0.3]) ;
    xlim([-N*dt 0]) ;

    for iisig = 1:nsigL
        subplot(nsigL+1,1,iisig+1) ;
        plot(tp,C*stcL(:,iisig)/norm(C*stcL(:,iisig)),clr{iisig},'LineWidth',3) ; hold on
        plot(tp,C*vdCL(:,isigL(iisig))/norm(C*vdCL(:,isigL(iisig))),'--','Color',clr{iisig},'LineWidth',3) ; hold on
        plot(tp,stc(:,iisig),'Color',[0.7 0.7 0.7],'LineWidth',3) ; hold on
        legend(['stcL' num2str(iisig) ' orth to staL'],['stcL' num2str(iisig) ', proj on staL = ' num2str(vdCL(:,isigL(iisig))'*staL,3)],['stc' num2str(iisig)]) ;
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
