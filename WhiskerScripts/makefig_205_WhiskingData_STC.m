% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_205 and generates plots that are
% shown in Figure 7A,B of the manuscript: 
% panel A - the spike triggered covariance features that were found to be
%           significant
% panel B - the eigenvalue of dC and the extent of o null eigenvalue
%           distribution used to determine which eigenvalues are
%           significant

% In addition we plot for every neuron:
% 1. stimulus distribution projected on the STA and STC features,
%    the spike conditioned stimulus distribution and the probability of a
%    spike (per time bin) given a stimulus (projected on a feature).
%    these are analogous to Figure 7C that shows results for whitened STC
% 2. Each STC feature is compared before and after projecting out the STA,
% 3. The entire stimulus set projected onto the STA-STC1 subspace, and the
%    stimuli that led to a spike
% 4. The 2D nonlinearity (analogous to Figure 7C that shows results for 
%    whitened STC)

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
    load([datadir 'VPM_cell_' Names{i} '_stcsig_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stc_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
    S = S/N ;
    iR = R == 1 ;
    
    z = zeros(T,2) ;
    z(:,1) = S*sta ;
    z(:,2) = S*stc(:,1) ;
    
    zx = max(abs(z(:))) ;
    
    nbnsplot = 100 ;
    Vzplot = linspace(-zx,zx,nbnsplot) ;
    [zplot1,zplot2] = meshgrid(Vzplot,Vzplot) ;
    zplot1 = reshape(zplot1,nbnsplot^2,1) ;
    zplot2 = reshape(zplot2,nbnsplot^2,1) ;
    Mzplot = [zplot1 zplot2] ;
    stc_surf = reshape(stc_model_rect_norm(Mzplot),nbnsplot,nbnsplot) ;
    
    
    ylm = 1.2*max(abs([edC(:) ; edCnull(:)]))*[-1 1] ;
    figure ;
    subplot(1,3,1) ;
    fill([0 0 N N 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
    plot(setdiff(1:N,isig),edC(setdiff(1:N,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
    plot(isig,edC(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
    plot(1:N,prctenull,'Color',0.9*ones(1,3),'LineWidth',1) ;
    ylim(ylm) ;
    xlim([0 N]) ;
    set(gca, 'Layer', 'top','FontSize',fntsz) ;
    title(lgnd{i},'FontSize',fntsz) ;
    ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
    xlabel('eigenvalue #','FontSize',fntsz) ;
    axis square ;
        
    subplot(1,3,2) ;
    plot(z(:,1),z(:,2),'.k') ; hold on ;
    plot(z(iR,1),z(iR,2),'.r') ;
    xlabel('S*sta','FontSize',fntsz) ;
    ylabel('S*stc_1','FontSize',fntsz) ;
    set(gca,'FontSize',fntsz) ;
    axis square xy ; box on ;
    xlim([-zx zx]) ;
    ylim([-zx zx]) ;
    
    subplot(1,3,3) ;
    imagesc(Vzplot,Vzplot,(stc_surf/dt),[0 180]) ;
    xlabel('S*sta','FontSize',fntsz) ;
    ylabel('S*stc_1','FontSize',fntsz) ;
    axis square xy ; box on ;
    colormap('cool') ;
    xlim([-zx zx]) ;
    ylim([-zx zx]) ;
    colorbar('eastoutside') ;
    set(gca,'FontSize',fntsz) ;
    
    figure ; 
    subplot(nsig+1,1,1) ;
    plot(tp,sta,'b','LineWidth',3) ; hold on 
    legend('sta') ;
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ; axis tight square ; box off
    xlabel('t (ms)','FontSize',fntsz) ;
    ylim([-0.3 0.3]) ;
    xlim([-N*dt 0]) ;
    
    for iisig = 1:nsig
        subplot(nsig+1,1,iisig+1) ;
    plot(tp,stc(:,iisig),clr{iisig},'LineWidth',3) ; hold on 
    plot(tp,vdC(:,isig(iisig)),'--','Color',clr{iisig},'LineWidth',3) ; hold on 
    legend(['stc' num2str(iisig) ' orth to sta'],['stc' num2str(iisig) ', proj on sta = ' num2str(vdC(:,isig(iisig))'*sta,3)]) ;
    set(gca,'FontSize',fntsz,'Yaxislocation','right') ; axis tight square ; box off
    xlabel('t (ms)','FontSize',fntsz) ;
    ylim([-0.3 0.3]) ;
    xlim([-N*dt 0]) ;

    end
    
                figure ;
                
                ctr_plot = bns_plot{1}(1:nbns_plot-1)+dbns_plot/2 ;
                ctr = bns{1}(1:nbns-1)+dbns/2 ;
                
                subplot(2,2,1) ;

                plot(ctr,Psf_a/dt,'b','LineWidth',2) ;
                axis tight square ; box off 
                ylim([0 0.25/dt]) ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                
                set(gca,'FontSize',fntsz) ;

                subplot(2,2,3) ;

                plot(ctr,Psf_c/dt,'b','LineWidth',2) ;
                axis tight square ; box off 
                ylim([0 0.25/dt]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('firing rate [Hz]','FontSize',fntsz) ;
                
                set(gca,'FontSize',fntsz) ;
                
                subplot(2,2,2) ;

                plot(ctr_plot,Pf_a_plot,'k','LineWidth',2) ; hold on ; 
                plot(ctr_plot,Pfs_a_plot,'r','LineWidth',2) ; hold on ;
                axis tight square ; box off 
                ylim([0 1]) ;
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('P(S*sta)','FontSize',fntsz) ;
                
                set(gca,'FontSize',fntsz) ;

                subplot(2,2,4) ;

                plot(ctr_plot,Pf_c_plot,'k','LineWidth',2) ; hold on ; 
                plot(ctr_plot,Pfs_c_plot,'r','LineWidth',2) ; hold on ;
                axis tight square ; box off 
                ylim([0 1]) ;
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('P(S*stc_1)','FontSize',fntsz) ;
                
                set(gca,'FontSize',fntsz) ;
end
