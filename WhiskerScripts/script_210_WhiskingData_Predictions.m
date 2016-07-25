% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script loads all the models that were fitted for each cell and runs
% them on the test portion of the stimulus to give a predicted spike train.

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

sr = 1000 ;        % (Hz) sampling rate
N = 150 ;          % stimulus dimensionality

nJK = 5 ;          % number of jackknives 

for i = 1:7
    for iJK = 1:nJK
        % Predicted spike train from STA model
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_sta_JK_' num2str(iJK) '.mat']) ;
        
        St = St/N ;
        Rt_sta = sta_model_rect_norm(St*sta) ;
        
        clearvars -except datadir Rt i Names N Rt_sta nJK iJK
        
        % Predicted spike train from Whitened-STA model
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_wsta_logl_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_wsta_L' num2str(Lopt) '_JK_' num2str(iJK) '.mat']) ;
        St = St/N ;
        StL = St*Cpinv ;
        
        Rt_staL = sta_model_rect_normL(StL*staL) ;
        
        clearvars -except datadir Rt i Names N Rt_sta Rt_staL nJK iJK
        
        % Predicted spike train from Whitened-STC model
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_wstc_logl_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_wstc_L' num2str(Lopt) '_JK_' num2str(iJK) '.mat']) ;
        
        St = St/N ;
        StL = St*Cpinv ;
        
        ztL = zeros(Tt,2) ;
        ztL(:,1) = StL*staL ;
        ztL(:,2) = StL*stcL(:,1) ;
        
        Rt_stcL = stc_model_rect_normL(ztL) ;
        Rt_stcL(isnan(Rt_stcL)) = mean(Rt_stcL(~isnan(Rt_stcL))) ;
        clearvars -except datadir Rt i Names N Rt_sta Rt_staL Rt_stcL nJK iJK
        
        % Predicted spike train from STC model
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_stc_JK_' num2str(iJK) '.mat']) ;
        
        St = St/N ;
        
        zt = zeros(Tt,2) ;
        zt(:,1) = St*sta ;
        zt(:,2) = St*stc(:,1) ;
        
        Rt_stc = stc_model_rect_norm(zt) ;
        Rt_stc(isnan(Rt_stc)) = mean(Rt_stc(~isnan(Rt_stc))) ;
        
        clearvars -except datadir Rt i Names N Rt_sta Rt_staL Rt_stc Rt_stcL nJK iJK
        
        % Predicted spike train from MNE model
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_mne_JK_' num2str(iJK) '.mat']) ;
        
        St = St - repmat(mean(St),[Tt,1]);
        St = St./repmat(std(St),[Tt,1]);
        
        SHt  = St*H' ;
        SJSt = sum(St.*(St*J),2) ;
        
        Rt_mne = 1./(1+exp(SJSt+SHt+A)) ;
        
        clearvars -except datadir St Spht Rt i Names N Rt_sta Rt_staL Rt_stc Rt_mne Rt_stcL nJK iJK
        
        % Predicted spike train from GLM (pre-computed in script 208)
        load([datadir 'VPM_cell_' Names{i} '_glm_JK_' num2str(iJK) '.mat']) ;
        clearvars -except datadir St Spht Rt i Names N Rt_sta Rt_staL Rt_stc Rt_mne Rt_glm Rt_stcL nJK iJK
        
        % Predicted spike train from tuning-curve model
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        load([datadir 'VPM_cell_' Names{i} '_tune_JK_' num2str(iJK) '.mat']) ;
        Rt_tun = tun_model_rect_norm(Spht) ;
        
        save([datadir 'VPM_cell_' Names{i} '_PredictedSpikeTrains_JK_' num2str(iJK) '.mat'],'Rt','Rt_sta','Rt_staL','Rt_stc','Rt_stcL','Rt_glm','Rt_mne','Rt_tun') ;
    end
end