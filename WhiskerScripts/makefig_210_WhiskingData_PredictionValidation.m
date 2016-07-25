% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_210, 211 and generates the plots 
% shown in Figures 11 and 12 of the manuscript:
% measured spike trains, predicted firing rates for all fitted models, log
% likelihood of each model (and the null model), power spectrum of
% stimulus, measured spike train and predicted response, coherence and 
% phase, coherence at specific frequency

clear ;
c_sta   = [  0   0 128]/255 ;
c_staL  = [128 128 255]/255 ;
c_stc   = [  0 128   0]/255 ;
c_stcL  = [128 255 128]/255 ;
c_mne   = [128   0 255]/255 ;
c_glm   = [255 128   0]/255 ;
c_tun   = [128 128 128]/255 ;
white   = [255 255 255]/255 ;
c_stm   = [0 0 0] ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

fmax=20;

nJK = 5 ; 

nfreq = 1000 ;
freqx = linspace(0,fmax,nfreq) ;
exmp_freq = [0 0 6.7 0 0.5 6 0] ;

SRt_JK        = zeros(nfreq,nJK) ;
SRt_sta_JK    = zeros(nfreq,nJK) ;
SRt_staL_JK   = zeros(nfreq,nJK) ;
SRt_stc_JK    = zeros(nfreq,nJK) ;
SRt_stcL_JK   = zeros(nfreq,nJK) ;
SRt_mne_JK    = zeros(nfreq,nJK) ;
SRt_glm_JK    = zeros(nfreq,nJK) ;
SRt_tun_JK    = zeros(nfreq,nJK) ;
SSt_JK        = zeros(nfreq,nJK) ;

Coh_sta_JK    = zeros(nfreq,nJK) ;
Coh_staL_JK   = zeros(nfreq,nJK) ;
Coh_stc_JK    = zeros(nfreq,nJK) ;
Coh_stcL_JK   = zeros(nfreq,nJK) ;
Coh_mne_JK    = zeros(nfreq,nJK) ;
Coh_glm_JK    = zeros(nfreq,nJK) ;
Coh_tun_JK    = zeros(nfreq,nJK) ;

phi_sta_JK    = zeros(nfreq,nJK) ;
phi_staL_JK   = zeros(nfreq,nJK) ;
phi_stc_JK    = zeros(nfreq,nJK) ;
phi_stcL_JK   = zeros(nfreq,nJK) ;
phi_mne_JK    = zeros(nfreq,nJK) ;
phi_glm_JK    = zeros(nfreq,nJK) ;
phi_tun_JK    = zeros(nfreq,nJK) ;


for i =  1:7
    for iJK = 1:nJK
        load([datadir 'VPM_cell_' Names{i} '_validation_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_loglikelihood.mat']) ;
        
        SRt_JK(:,iJK)        = interp1(freq,SRt,freqx) ;
        SRt_sta_JK(:,iJK)    = interp1(freq,SRt_sta,freqx) ;
        SRt_staL_JK(:,iJK)   = interp1(freq,SRt_staL,freqx) ;
        SRt_stc_JK(:,iJK)    = interp1(freq,SRt_stc,freqx) ;
        SRt_stcL_JK(:,iJK)   = interp1(freq,SRt_stcL,freqx) ;
        SRt_mne_JK(:,iJK)    = interp1(freq,SRt_mne,freqx) ;
        SRt_glm_JK(:,iJK)    = interp1(freq,SRt_glm,freqx) ;
        SRt_tun_JK(:,iJK)    = interp1(freq,SRt_tun,freqx) ;
        SSt_JK(:,iJK)        = interp1(freq,SSt,freqx) ;

        Coh_sta_JK(:,iJK)    = interp1(freq,Coh_sta,freqx) ;
        Coh_staL_JK(:,iJK)   = interp1(freq,Coh_staL,freqx) ;
        Coh_stc_JK(:,iJK)    = interp1(freq,Coh_stc,freqx) ;
        Coh_stcL_JK(:,iJK)   = interp1(freq,Coh_stcL,freqx) ;
        Coh_mne_JK(:,iJK)    = interp1(freq,Coh_mne,freqx) ;
        Coh_glm_JK(:,iJK)    = interp1(freq,Coh_glm,freqx) ;
        Coh_tun_JK(:,iJK)    = interp1(freq,Coh_tun,freqx) ;
        
        phi_sta_JK(:,iJK)    = interp1(freq,phi_sta_s,freqx) ;
        phi_staL_JK(:,iJK)   = interp1(freq,phi_staL_s,freqx) ;
        phi_stc_JK(:,iJK)    = interp1(freq,phi_stc_s,freqx) ;
        phi_stcL_JK(:,iJK)   = interp1(freq,phi_stcL_s,freqx) ;
        phi_mne_JK(:,iJK)    = interp1(freq,phi_mne_s,freqx) ;
        phi_glm_JK(:,iJK)    = interp1(freq,phi_glm_s,freqx) ;
        phi_tun_JK(:,iJK)    = interp1(freq,phi_tun_s,freqx) ;
    
    end
    
    m_SRt      = mean(SRt_JK,2) ;
    m_SRt_sta  = mean(SRt_sta_JK,2) ;
    m_SRt_staL = mean(SRt_staL_JK,2) ;
    m_SRt_stc  = mean(SRt_stc_JK,2) ;
    m_SRt_stcL = mean(SRt_stcL_JK,2) ;
    m_SRt_mne  = mean(SRt_mne_JK,2) ;
    m_SRt_glm  = mean(SRt_glm_JK,2) ;
    m_SRt_tun  = mean(SRt_tun_JK,2) ;
    m_SSt      = mean(SSt_JK,2) ;

    s_SRt      = std(SRt_JK,[],2)/sqrt(nJK) ;
    s_SRt_sta  = std(SRt_sta_JK,[],2)/sqrt(nJK) ;
    s_SRt_staL = std(SRt_staL_JK,[],2)/sqrt(nJK) ;
    s_SRt_stc  = std(SRt_stc_JK,[],2)/sqrt(nJK) ;
    s_SRt_stcL = std(SRt_stcL_JK,[],2)/sqrt(nJK) ;
    s_SRt_mne  = std(SRt_mne_JK,[],2)/sqrt(nJK) ;
    s_SRt_glm  = std(SRt_glm_JK,[],2)/sqrt(nJK) ;
    s_SRt_tun  = std(SRt_tun_JK,[],2)/sqrt(nJK) ;
    s_SSt      = std(SSt_JK,[],2)/sqrt(nJK) ;
    
    m_Coh_sta  = mean(Coh_sta_JK,2) ;
    m_Coh_staL = mean(Coh_staL_JK,2) ;
    m_Coh_stc  = mean(Coh_stc_JK,2) ;
    m_Coh_stcL = mean(Coh_stcL_JK,2) ;
    m_Coh_mne  = mean(Coh_mne_JK,2) ;
    m_Coh_glm  = mean(Coh_glm_JK,2) ;
    m_Coh_tun  = mean(Coh_tun_JK,2) ;
  
    s_Coh_sta  = std(Coh_sta_JK,[],2)/sqrt(nJK) ;
    s_Coh_staL = std(Coh_staL_JK,[],2)/sqrt(nJK) ;
    s_Coh_stc  = std(Coh_stc_JK,[],2)/sqrt(nJK) ;
    s_Coh_stcL = std(Coh_stcL_JK,[],2)/sqrt(nJK) ;
    s_Coh_mne  = std(Coh_mne_JK,[],2)/sqrt(nJK) ;
    s_Coh_glm  = std(Coh_glm_JK,[],2)/sqrt(nJK) ;
    s_Coh_tun  = std(Coh_tun_JK,[],2)/sqrt(nJK) ;
 
    m_phi_sta  = mean(phi_sta_JK,2) ;
    m_phi_staL = mean(phi_staL_JK,2) ;
    m_phi_stc  = mean(phi_stc_JK,2) ;
    m_phi_stcL = mean(phi_stcL_JK,2) ;
    m_phi_mne  = mean(phi_mne_JK,2) ;
    m_phi_glm  = mean(phi_glm_JK,2) ;
    m_phi_tun  = mean(phi_tun_JK,2) ;
  
    s_phi_sta  = std(phi_sta_JK,[],2)/sqrt(nJK) ;
    s_phi_staL = std(phi_staL_JK,[],2)/sqrt(nJK) ;
    s_phi_stc  = std(phi_stc_JK,[],2)/sqrt(nJK) ;
    s_phi_stcL = std(phi_stcL_JK,[],2)/sqrt(nJK) ;
    s_phi_mne  = std(phi_mne_JK,[],2)/sqrt(nJK) ;
    s_phi_glm  = std(phi_glm_JK,[],2)/sqrt(nJK) ;
    s_phi_tun  = std(phi_tun_JK,[],2)/sqrt(nJK) ;
 
    [~,ifreq] = min(abs(freqx-exmp_freq(i))) ; 
    
    figure ;

    subplot(3,4,1:3) ;
    plot(freqx,log(m_SRt),'Color',[200 200 200]/255) ; hold on ;
    plot(freqx,log(m_SRt_sta),'Color',c_sta,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SRt_staL),'Color',c_staL,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SRt_stc),'Color',c_stc,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SRt_stcL),'Color',c_stcL,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SRt_mne),'Color',c_mne,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SRt_glm),'Color',c_glm,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SRt_tun),'Color',c_tun,'LineWidth',2) ; hold on ;
    plot(freqx,log(m_SSt),'Color',c_stm,'LineWidth',2) ; hold on ;
    plot(freqx(ifreq),-2,'*k','MarkerSize',15) ; hold on ;
    ylabel('Log(Spectral Power)','FontSize',10) ;
    title(['cell ' lgnd{i} ': power spectrum and coherence']) ;
    
    subplot(3,4,5:7);
    plot(freqx,(m_Coh_sta),'Color',c_sta,'LineWidth',2) ; hold on ;
    plot(freqx,(m_Coh_staL),'Color',c_staL,'LineWidth',2) ; hold on ;
    plot(freqx,(m_Coh_stc),'Color',c_stc,'LineWidth',2) ; hold on ;
    plot(freqx,(m_Coh_stcL),'Color',c_stcL,'LineWidth',2) ; hold on ;
    plot(freqx,(m_Coh_mne),'Color',c_mne,'LineWidth',2) ; hold on ;
    plot(freqx,(m_Coh_glm),'Color',c_glm,'LineWidth',2) ; hold on ;
    plot(freqx,(m_Coh_tun),'Color',c_tun,'LineWidth',2) ; hold on ;
    plot(freqx(ifreq),1,'*k','MarkerSize',15) ; hold on ;
    plot([0 fmax],confC*ones(1,2),'Color',c_sta,'LineWidth',2) ; hold on ; % plot confidence interval
    ylabel('Magnitude of Coherence - |C|','FontSize',10) ;
    ylim([0 1]) ;

    subplot(3,4,9:11); 
    
    plot(freqx,m_phi_sta,'Color',c_sta,'LineWidth',2) ; hold on ;
    plot(freqx,m_phi_staL,'Color',c_staL,'LineWidth',2) ; hold on ;
    plot(freqx,m_phi_stc,'Color',c_stc,'LineWidth',2) ; hold on ;
    plot(freqx,m_phi_stcL,'Color',c_stcL,'LineWidth',2) ; hold on ;
    plot(freqx,m_phi_mne,'Color',c_mne,'LineWidth',2) ; hold on ;
    plot(freqx,m_phi_glm,'Color',c_glm,'LineWidth',2) ; hold on ;
    plot(freqx,m_phi_tun,'Color',c_tun,'LineWidth',2) ; hold on ;
    
    ylabel('Phase of Coherence - arg{C}','FontSize',10) ;
    set(gca,'ytick',-pi:pi/2:pi,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}) ;
    ylim([-pi pi]) ;
    xlim([0 fmax]) ;
    
    m_LL = mean(LL(i,:,:),3) ;
    s_LL = std(LL(i,:,:),[],3)/sqrt(nJK) ;
    
    subplot(3,4,12) ;
    bar(1,m_LL(1),'FaceColor',c_sta) ; hold on ;
    bar(2,m_LL(2),'FaceColor',c_staL) ; hold on ;
    bar(3,m_LL(3),'FaceColor',c_stc) ; hold on ;
    bar(4,m_LL(4),'FaceColor',c_stcL) ; hold on ;
    bar(5,m_LL(5),'FaceColor',c_mne) ; hold on ;
    bar(6,m_LL(6),'FaceColor',c_glm) ; hold on ;
    bar(7,m_LL(7),'FaceColor',c_tun) ; hold on ;
    bar(8,m_LL(8),'FaceColor','w') ; hold on ;
    ylabel('log likelihood') ;
    set(gca,'xtick',1:8,'xticklabel',{'sta','pwsta','stc','pwstc','mne','glm','tuning','null'},'xticklabelrotation',90) ;
    xlim([0.5 8.5]) ;
    for imod = 1:8
        plot([imod imod],m_LL(imod)+s_LL(imod)*[-1 1],'k','LineWidth',2) ; hold on ;
    end

    % plot(1:8,squeeze(LL(i,:,:)),'--ok','MarkerSize',3,'MarkerFaceColor','k') ;
    
    subplot(3,4,8) ;
    bar(1,m_Coh_sta(ifreq),'FaceColor',c_sta) ; hold on ;
    bar(2,m_Coh_staL(ifreq),'FaceColor',c_staL) ; hold on ;
    bar(3,m_Coh_stc(ifreq),'FaceColor',c_stc) ; hold on ;
    bar(4,m_Coh_stcL(ifreq),'FaceColor',c_stcL) ; hold on ;
    bar(5,m_Coh_mne(ifreq),'FaceColor',c_mne) ; hold on ;
    bar(6,m_Coh_glm(ifreq),'FaceColor',c_glm) ; hold on ;
    bar(7,m_Coh_tun(ifreq),'FaceColor',c_tun) ; hold on ;
    
    plot([1 1],m_Coh_sta(ifreq)+s_Coh_sta(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    plot([2 2],m_Coh_staL(ifreq)+s_Coh_staL(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    plot([3 3],m_Coh_stc(ifreq)+s_Coh_stc(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    plot([4 4],m_Coh_stcL(ifreq)+s_Coh_stcL(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    plot([5 5],m_Coh_mne(ifreq)+s_Coh_mne(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    plot([6 6],m_Coh_glm(ifreq)+s_Coh_glm(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    plot([7 7],m_Coh_tun(ifreq)+s_Coh_tun(ifreq)*[-1 1],'k','LineWidth',2) ; hold on ;
    ylabel(['coherence at f = ' num2str(exmp_freq(i))]) ;
    set(gca,'xtick',1:7,'xticklabel',{'sta','pwsta','stc','pwstc','mne','glm','tuning'},'xticklabelrotation',90) ;
    xlim([0.5 8.5]) ;
    ylim([0 1]) ;
    
    load([datadir 'VPM_cell_' Names{i} '_PredictedSpikeTrains_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        
    figure ; 
    stem(t,Rt     ,'Color',[200 200 200]/255,'LineWidth',1,'Marker','None') ; hold on ;
    plot(t,Rt_sta ,'Color',c_sta,'LineWidth',2) ; hold on ;
    plot(t,Rt_staL,'Color',c_staL,'LineWidth',2) ; hold on ;
    plot(t,Rt_stc ,'Color',c_stc,'LineWidth',2) ; hold on ;
    plot(t,Rt_stcL,'Color',c_stcL,'LineWidth',2) ; hold on ;
    plot(t,Rt_mne ,'Color',c_mne,'LineWidth',2) ; hold on ;
    plot(t,Rt_glm ,'Color',c_glm,'LineWidth',2) ; hold on ;
    plot(t,Rt_tun ,'Color',c_tun,'LineWidth',2) ; hold on ;

end