% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads the output of script_108, 109 and generates the plots 
% shown in Figure 10 of the manuscript:
% panel A - measured spike train
% panel B - predicted firing rates for STA, STC and GLM models
% panel C - log likelihood of each model (and the null model)
% panel D - power spectrum of stimulus, measured spike train and predicted
%           response
% panel E - coherence and phase
% panel F - coherence at low frequencies

clear ;

icell = 3 ;
iL = 3 ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

c_sta   = [  0   0 128]/255 ;
c_stc   = [  0 128   0]/255 ;
c_glm   = [255 128   0]/255 ;
c_nll   = [255 255 255]/255 ;
c_stm   = [0 0 0] ;

stim_length = {'short2','short3','long'} ;

load([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_validation_JK_1.mat']) ;
load([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_loglikelihood.mat']) ;

nF = length(freq) ;
nJK = 5 ;

Coh_sta_JK = zeros(nF,nJK) ;
Coh_stc_JK = zeros(nF,nJK) ;
Coh_glm_JK = zeros(nF,nJK) ;

phi_sta_JK = zeros(nF,nJK) ;
phi_stc_JK = zeros(nF,nJK) ;
phi_glm_JK = zeros(nF,nJK) ;

SRt_JK     = zeros(nF,nJK) ;
SRt_sta_JK = zeros(nF,nJK) ;
SRt_stc_JK = zeros(nF,nJK) ;
SRt_glm_JK = zeros(nF,nJK) ;

LL_JK = LL ;

for iJK = 1:nJK 
    load([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_validation_JK_' num2str(iJK) '.mat']) ;
    
    Coh_sta_JK(:,iJK) = mean(Coh_sta,2) ;
    Coh_stc_JK(:,iJK) = mean(Coh_stc,2) ;
    Coh_glm_JK(:,iJK) = mean(Coh_glm,2) ;
    
    phi_sta_JK(:,iJK) = mean(phi_sta_s,2) ;
    phi_stc_JK(:,iJK) = mean(phi_stc_s,2) ;
    phi_glm_JK(:,iJK) = mean(phi_glm_s,2) ;
    
    SRt_JK    (:,iJK) = SRt     ;
    SRt_sta_JK(:,iJK) = SRt_sta ;
    SRt_stc_JK(:,iJK) = SRt_stc ;
    SRt_glm_JK(:,iJK) = SRt_glm ;
end

m_Coh_sta = mean(Coh_sta_JK,2) ;
m_Coh_stc = mean(Coh_stc_JK,2) ;
m_Coh_glm = mean(Coh_glm_JK,2) ;

m_phi_sta = mean(phi_sta_JK,2) ;
m_phi_stc = mean(phi_stc_JK,2) ;
m_phi_glm = mean(phi_glm_JK,2) ;

m_SRt     = mean(SRt_JK    ,2) ;
m_SRt_sta = mean(SRt_sta_JK,2) ;
m_SRt_stc = mean(SRt_stc_JK,2) ;
m_SRt_glm = mean(SRt_glm_JK,2) ;

s_Coh_sta = std(Coh_sta_JK,[],2)/sqrt(nJK) ;
s_Coh_stc = std(Coh_stc_JK,[],2)/sqrt(nJK) ;
s_Coh_glm = std(Coh_glm_JK,[],2)/sqrt(nJK) ;

s_phi_sta = std(phi_sta_JK,[],2)/sqrt(nJK) ;
s_phi_stc = std(phi_stc_JK,[],2)/sqrt(nJK) ;
s_phi_glm = std(phi_glm_JK,[],2)/sqrt(nJK) ;

s_SRt     = std(SRt_JK    ,[],2)/sqrt(nJK) ;
s_SRt_sta = std(SRt_sta_JK,[],2)/sqrt(nJK) ;
s_SRt_stc = std(SRt_stc_JK,[],2)/sqrt(nJK) ;
s_SRt_glm = std(SRt_glm_JK,[],2)/sqrt(nJK) ;

m_LL = mean(LL_JK,2) ;
s_LL = std(LL_JK,[],2)/sqrt(nJK) ;


figure ;

subplot(3,1,1) ;

plot(freq,log(m_SRt),'Color',[200 200 200]/255) ; hold on ;
plot(freq,log(m_SRt_sta),'Color',c_sta,'LineWidth',2) ; hold on ;
plot(freq,log(m_SRt_stc),'Color',c_stc,'LineWidth',2) ; hold on ;
plot(freq,log(m_SRt_glm),'Color',c_glm,'LineWidth',2) ; hold on ;
plot(freq,log(SSt),'k','LineWidth',2) ; hold on ;
plot([0.5 5],-4.5*[1 1],'-','LineWidth',5,'Color',[0.7 0.7 0.7]) ;
ylabel('Log(Spectral Power)','FontSize',10) ;
xlim([0 freq(nF)]) ;

subplot(3,1,2);
plot(freq,m_Coh_sta,'Color',c_sta,'LineWidth',2) ; hold on ;
plot(freq,m_Coh_stc,'Color',c_stc,'LineWidth',2) ; hold on ;
plot(freq,m_Coh_glm,'Color',c_glm,'LineWidth',2) ; hold on ;
plot([min(freq) max(freq)],confC*ones(1,2),'--k','LineWidth',2) ; hold on ; % plot confidence interval
plot([0.5 5],0.6*[1 1],'-','LineWidth',5,'Color',[0.7 0.7 0.7]) ;
ylabel('Magnitude of Coherence - |C|','FontSize',10) ;
ylim([0 1]) ;
xlim([0 freq(nF)]) ;

subplot(3,1,3);

plot(freq,m_phi_sta,'Color',c_sta,'LineWidth',2) ; hold on ;
plot(freq,m_phi_stc,'Color',c_stc,'LineWidth',2) ; hold on ;
plot(freq,m_phi_glm,'Color',c_glm,'LineWidth',2) ; hold on ;

ylabel('Phase of Coherence - arg{C}','FontSize',10) ;
set(gca,'ytick',-pi:pi/2:pi,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}) ;
ylim([-pi pi]) ;
xlim([0 freq(nF)]) ;

% NW = 60 
% Tt = 823 sec
bw_samp = 2*60/914 ; 
bw_mean = 4.5 ;
ifreq = intersect(find(freq>=0.5),find(freq<=5)) ; 

Coh_sta_JK_bnd = Coh_sta_JK(ifreq,:) ; 
Coh_stc_JK_bnd = Coh_stc_JK(ifreq,:) ; 
Coh_glm_JK_bnd = Coh_glm_JK(ifreq,:) ; 

m_Coh_sta_bnd = mean(Coh_sta_JK_bnd(:)) ; 
m_Coh_stc_bnd = mean(Coh_stc_JK_bnd(:)) ; 
m_Coh_glm_bnd = mean(Coh_glm_JK_bnd(:)) ; 

s_Coh_sta_bnd = std(Coh_sta_JK_bnd(:))/sqrt(nJK)/sqrt(bw_mean/bw_samp) ; 
s_Coh_stc_bnd = std(Coh_stc_JK_bnd(:))/sqrt(nJK)/sqrt(bw_mean/bw_samp) ; 
s_Coh_glm_bnd = std(Coh_glm_JK_bnd(:))/sqrt(nJK)/sqrt(bw_mean/bw_samp) ; 

figure ;
subplot(1,2,1) ;
bar(1,m_LL(1),'EdgeColor','None','FaceColor',c_sta) ; hold on ;
bar(2,m_LL(2),'EdgeColor','None','FaceColor',c_stc) ; hold on ;
bar(3,m_LL(3),'EdgeColor','None','FaceColor',c_glm) ; hold on ;
bar(4,m_LL(4),'EdgeColor','k','FaceColor',[1 1 1]) ; hold on ;

plot([1 1],m_LL(1) + s_LL(1)*[-1 1],'k','LineWidth',3) ; hold on ;
plot([2 2],m_LL(2) + s_LL(2)*[-1 1],'k','LineWidth',3) ; hold on ;
plot([3 3],m_LL(3) + s_LL(3)*[-1 1],'k','LineWidth',3) ; hold on ;
plot([4 4],m_LL(4) + s_LL(4)*[-1 1],'k','LineWidth',3) ; hold on ;

set(gca,'xtick',1:4,'xticklabel',{'sta','sta+stc','glm','constant'}) ;
ylabel('log likelihood') ;

subplot(1,2,2) ;

bar(1,m_Coh_sta_bnd,'EdgeColor','None','FaceColor',c_sta) ; hold on ;
bar(2,m_Coh_stc_bnd,'EdgeColor','None','FaceColor',c_stc) ; hold on ;
bar(3,m_Coh_glm_bnd,'EdgeColor','None','FaceColor',c_glm) ; hold on ;

plot([1 1],m_Coh_sta_bnd + s_Coh_sta_bnd*[-1 1],'k','LineWidth',3) ; hold on ;
plot([2 2],m_Coh_stc_bnd + s_Coh_stc_bnd*[-1 1],'k','LineWidth',3) ; hold on ;
plot([3 3],m_Coh_glm_bnd + s_Coh_glm_bnd*[-1 1],'k','LineWidth',3) ; hold on ;

set(gca,'xtick',1:3,'xticklabel',{'sta','sta+stc','glm'}) ;
ylabel(['average coherence at 0.5-5 Hz frequency band']) ;


load([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_Predictions_JK_' num2str(iJK) '.mat']) ;
load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
          
figure ;

stem(t,Rt     ,'Color',[200 200 200]/255,'Marker','None','LineWidth',2) ; hold on ;
plot(t,Rt_sta ,'Color',c_sta,'LineWidth',2) ; hold on ;
plot(t,Rt_stc ,'Color',c_stc,'LineWidth',2) ; hold on ;
plot(t,Rt_glm ,'Color',c_glm,'LineWidth',2) ; hold on ;

xlabel('t [s]') ;
ylabel('# of spikes') ;
legend('data','sta','stc','glm') ;
axis tight
