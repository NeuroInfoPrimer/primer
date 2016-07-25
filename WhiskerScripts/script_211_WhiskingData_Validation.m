% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the validation metrics of the predicted spike trains
% from the different analysis methods/models. 
% these metrics are the log-likelihood and the coherence
% coherence is computed using the chronux toolbox which should appear in
% the MATLAB search path. In addition to the coherence itself, this script
% computes the spectrum of the stimulus, the measured response and the
% simulated response as well as the coherence phase for each model (phase
% with respect to the recorded spike train).

% for details see Equations 48-50 and Box 7

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

if exist('coherencyc','file') == 0
    error('this script uses the Chronux toolbox. please add the toolbox to the search path and rerun') ;
end

nJK = 5 ; % number of jackknives

LL = zeros(7,8,nJK) ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;

sr = 1000 ;        % (Hz) sampling rate
ds = 2 ;           % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 

fmax=25 ; % [Hz] maximal frequency in spectrum computations

for i = 1:7
    
    for iJK = 1:nJK
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ;
        % the stimulus has one spatial dimension, so only one of the 
        % stimulus coordinates is used to compute the spectrum
        St = St(:,1) ; 
    
        load([datadir 'VPM_cell_' Names{i} '_PredictedSpikeTrains_JK_' num2str(iJK) '.mat']) ;
        t = dt*(1:Tt) ; % just the time base

        % implementation of Equation 48
        LL(i,1,iJK) = mean(Rt.*log(Rt_sta)-Rt_sta) ;
        LL(i,2,iJK) = mean(Rt.*log(Rt_staL)-Rt_staL) ;
        LL(i,3,iJK) = mean(Rt.*log(Rt_stc)-Rt_stc) ;
        LL(i,4,iJK) = mean(Rt.*log(Rt_stcL)-Rt_stcL) ;
        LL(i,5,iJK) = mean(Rt.*log(Rt_mne)-Rt_mne) ;
        LL(i,6,iJK) = mean(Rt.*log(Rt_glm)-Rt_glm) ;
        LL(i,7,iJK) = mean(Rt.*log(Rt_tun)-Rt_tun) ;
        
        % implementation of Equation 49
        LL(i,8,iJK) = mean(Rt.*log(mean(Rt))-mean(Rt)) ;
        

        NW=12 ;                     % time-bandwidth product (for spectrum calculation)
        params.tapers = [NW 2*NW-1];
        params.Fs = 1/dt ;          % sampling frequency
        params.fpass=[0 fmax];      % frequency range
        params.pad = 4;             % 4-times padding, a rule of thumb
        params.trialave = 0;        % one trace
    
        % using chronux routines to compute spectrum
        [SRt,     ~] = mtspectrumc(Rt-mean(Rt),params); 
        [SRt_sta, ~] = mtspectrumc(Rt_sta-mean(Rt_sta),params);
        [SRt_staL,~] = mtspectrumc(Rt_staL-mean(Rt_staL),params);
        [SRt_stc, ~] = mtspectrumc(Rt_stc-mean(Rt_stc),params);
        [SRt_stcL,~] = mtspectrumc(Rt_stcL-mean(Rt_stcL),params);
        [SRt_mne, ~] = mtspectrumc(Rt_mne-mean(Rt_mne),params);
        [SRt_glm, ~] = mtspectrumc(Rt_glm-mean(Rt_glm),params);
        [SRt_tun, ~] = mtspectrumc(Rt_tun-mean(Rt_tun),params);
        [SSt    , ~] = mtspectrumc(St(:,1)-mean(St(:,1)),params);
    
        SRt      = SRt/norm(SRt) ;
        SRt_sta  = SRt_sta/norm(SRt_sta) ;
        SRt_staL = SRt_staL/norm(SRt_staL) ;
        SRt_stc  = SRt_stc/norm(SRt_stc) ;
        SRt_stcL = SRt_stcL/norm(SRt_stcL) ;
        SRt_mne  = SRt_mne/norm(SRt_mne) ;
        SRt_glm  = SRt_glm/norm(SRt_glm) ;
        SRt_tun  = SRt_tun/norm(SRt_tun) ;
        SSt      = SSt/norm(SSt) ;
    
        NW=25 ;                     % time-bandwidth product (for coherence calculation)
        params.tapers = [NW 2*NW-1] ;
        params.Fs = 1/dt ;
        params.fpass=[0 fmax] ;
        params.pad = 4 ;
        params.err= [2 0.05] ;      % confidence interval
        params.trialave = 0 ;

        % using chronux routines to compute coherence, phase of predicted firing rates and measured response
        [Coh_sta,phi_sta,~,~,~,freq,confC,~,~] = coherencyc(Rt-mean(Rt),Rt_sta-mean(Rt_sta),params) ;
        [Coh_staL,phi_staL,~,~,~,~,~,~,~]      = coherencyc(Rt-mean(Rt),Rt_staL-mean(Rt_staL),params) ;
        [Coh_stc,phi_stc,~,~,~,~,~,~,~]        = coherencyc(Rt-mean(Rt),Rt_stc-mean(Rt_stc),params) ;
        [Coh_stcL,phi_stcL,~,~,~,~,~,~,~]      = coherencyc(Rt-mean(Rt),Rt_stcL-mean(Rt_stcL),params) ;
        [Coh_mne,phi_mne,~,~,~,~,~,~,~]        = coherencyc(Rt-mean(Rt),Rt_mne-mean(Rt_mne),params) ;
        [Coh_glm,phi_glm,~,~,~,~,~,~,~]        = coherencyc(Rt-mean(Rt),Rt_glm-mean(Rt_glm),params) ;
        [Coh_tun,phi_tun,~,~,~,~,~,~,~]        = coherencyc(Rt-mean(Rt),Rt_tun-mean(Rt_tun),params) ; 
        
        % phase is computed modulo 2pi. in frequencies where the coherence is
        % not found to be statistically significant, the phase is set to NaN
        phi_sta_s = mod(phi_sta,2*pi) ; 
        phi_sta_s(phi_sta_s>pi) = phi_sta_s(phi_sta_s>pi) - 2*pi ; 
        phi_sta_s(Coh_sta<confC) = NaN ;
    
        phi_staL_s = mod(phi_staL,2*pi) ; 
        phi_staL_s(phi_staL_s>pi) = phi_staL_s(phi_staL_s>pi) - 2*pi ; 
        phi_staL_s(Coh_staL<confC) = NaN ;
    
        phi_stc_s = mod(phi_stc,2*pi) ; 
        phi_stc_s(phi_stc_s>pi) = phi_stc_s(phi_stc_s>pi) - 2*pi ; 
        phi_stc_s(Coh_stc<confC) = NaN ;
    
        phi_stcL_s = mod(phi_stcL,2*pi) ; 
        phi_stcL_s(phi_stcL_s>pi) = phi_stcL_s(phi_stcL_s>pi) - 2*pi ; 
        phi_stcL_s(Coh_stcL<confC) = NaN ;
        
        phi_mne_s = mod(phi_mne,2*pi) ;
        phi_mne_s(phi_mne_s>pi) = phi_mne_s(phi_mne_s>pi) - 2*pi ;
        phi_mne_s(Coh_mne<confC) = NaN ;
        
        phi_glm_s = mod(phi_glm,2*pi) ;
        phi_glm_s(phi_glm_s>pi) = phi_glm_s(phi_glm_s>pi) - 2*pi ;
        phi_glm_s(Coh_glm<confC) = NaN ;
        
        phi_tun_s = mod(phi_tun,2*pi) ;
        phi_tun_s(phi_tun_s>pi) = phi_tun_s(phi_tun_s>pi) - 2*pi ;
        phi_tun_s(Coh_tun<confC) = NaN ;
        save([datadir 'VPM_cell_' Names{i} '_validation_JK_' num2str(iJK) '.mat'],'t','freq','confC','SRt','SRt_sta','SRt_staL','SRt_stc','SRt_stcL','SRt_mne','SRt_glm','SRt_tun','SSt',...
                                                                        'Coh_sta','Coh_staL','Coh_stc','Coh_stcL','Coh_mne','Coh_glm','Coh_tun',...
                                                                        'phi_sta_s','phi_staL_s','phi_stc_s','phi_stcL_s','phi_mne_s','phi_glm_s','phi_tun_s') ;

    end
    save([datadir 'VPM_cell_' Names{i} '_loglikelihood.mat'],'LL') ;
end