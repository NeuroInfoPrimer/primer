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
datadir = 'RetinaData/' ;

if exist('coherencyc','file') == 0
    error('this script uses the Chronux toolbox. please add the toolbox to the search path and rerun') ;
end


icell = 3 ;
iL = 3 ;

dt = 1/30 ;       % delta t of stimulus
sr = 1/dt ;       % (Hz) sampling rate

fmax = 15 ; % [Hz] maximal frequency in spectrum computations

stim_length = {'short2','short3','long'} ;

nj = 4 ;    % number of chunks the responses (model and experiment) are
            % are split into. the spectral calculations may become heavy
            % for a long recording if it is not split

njs = 40 ;  % number of stimulus pixels for which the spectrum is computed
            % for this dataset the stimulus is a checkerboard stimulus with
            % no spatial or temporal correlations, so there is no need to
            % compute the spectrum for all stimulus dimensions. for other
            % stimuli it may be beneficial to do so. 

nJK = 5 ;   % number of jackknives

LL = zeros(4,nJK) ; 
            % variable that holds the log-likelihood values for 3 models
            % (STA, STC, GLM) (Equation 48) and the LL of the null model
            % (Equation 49), for each jack-knife (i.e., every 20% chunk of 
            % portion of the experiment the experiment that predictions
            % were generated for).

for iJK =  1:nJK 
    load([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_Predictions_JK_' num2str(iJK) '.mat']) ;
    load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;

    Tt = length(Rt) ;    % number of time bins for which a response was predicted
    Ttj = floor(Tt/nj) ; % number of time bins in each chunk of response for which the spectrum will be computed

    t = dt*(1:Tt) ;    % time in units of seconds    
    tj = dt*(1:Ttj) ;  % time in units of seconds for the chunk of response for which the spectrum will be computed
    
    Rtj = zeros(Ttj,nj) ;     % chunk of measured response
    Rtj_sta = zeros(Ttj,nj) ; % chunk of predicted response (STA model)
    Rtj_stc = zeros(Ttj,nj) ; % chunk of predicted response (STC model)
    Rtj_glm = zeros(Ttj,nj) ; % chunk of predicted response (GLM model)
    Stj = zeros(Ttj,nj*njs) ; % chunk of stimulus

    % implementation of Equation 48
    LL(1,iJK) = mean(Rt.*log(Rt_sta)-Rt_sta) ;
    LL(2,iJK) = mean(Rt.*log(Rt_stc)-Rt_stc) ;
    LL(3,iJK) = mean(Rt.*log(Rt_glm)-Rt_glm) ; 
    
    % implementation of Equation 49
    LL(4,iJK) = mean(Rt.*log(mean(Rt))-mean(Rt)) ;
    
    % spliting the response (model and experiment) and stimulus to nj and 
    % nj*njs parts respectively so that computation of spectrum is faster
    for j = 1:nj
        Rtj(:,j) = Rt((j-1)*Ttj+(1:Ttj))-mean(Rt((j-1)*Ttj+(1:Ttj))) ;
        Rtj_sta(:,j) = Rt_sta((j-1)*Ttj+(1:Ttj))-mean(Rt_sta((j-1)*Ttj+(1:Ttj))) ;
        Rtj_stc(:,j) = Rt_stc((j-1)*Ttj+(1:Ttj))-mean(Rt_stc((j-1)*Ttj+(1:Ttj))) ;
        Rtj_glm(:,j) = Rt_glm((j-1)*Ttj+(1:Ttj))-mean(Rt_glm((j-1)*Ttj+(1:Ttj))) ;
        for js = 1:njs
            Stj(:,(j-1)*njs+js) = St((j-1)*Ttj+(1:Ttj),js) ;
        end
    end

    NW = 80 ;                   % time-bandwidth product (for spectrum calculation)
    params.tapers = [NW 2*NW-1];  
    params.Fs = 1/dt ;          % sampling frequency
    params.fpass=[0 fmax] ;     % frequency range
    params.pad = 4 ;            % 4-times padding, a rule of thumb
    params.trialave = 1 ;        


    % using chronux routines to compute spectrum
    [SRt, ~] = mtspectrumc(Rtj,params); 
    [SRt_sta, ~] = mtspectrumc(Rtj_sta,params);
    [SRt_stc, ~] = mtspectrumc(Rtj_stc,params);
    [SRt_glm, ~] = mtspectrumc(Rtj_glm,params);
    SSt = zeros(length(SRt),njs) ;
    for js = 1:njs
        [SSt(:,js), ~] = mtspectrumc(Stj(:,(js-1)*nj+(1:nj)),params) ;
    end
    SSt = mean(SSt,2) ;
    
    % normalizing spectrum
    SRt = SRt/norm(SRt) ;
    SRt_sta = SRt_sta/norm(SRt_sta) ;
    SRt_stc = SRt_stc/norm(SRt_stc) ;
    SRt_glm = SRt_glm/norm(SRt_glm) ;
    SSt = SSt/norm(SSt) ;

    NW=60 ;                     % time-bandwidth product (for coherence calculation)
    params.tapers = [NW 2*NW-1] ;
    params.Fs = 1/dt ; 
    params.fpass=[0 fmax] ; 
    params.pad = 4 ; 
    params.err=[2 0.05] ;       % confidence interval 
    params.trialave = 1 ;
    
    % using chronux routines to compute coherence, phase of predicted firing rates and measured response
    [Coh_sta,phi_sta,~,~,~,freq,confC,~,~] = coherencyc(Rtj,Rtj_sta,params) ;
    [Coh_stc,phi_stc,~,~,~,~,~,~,~] = coherencyc(Rtj,Rtj_stc,params) ;
    [Coh_glm,phi_glm,~,~,~,~,~,~,~] = coherencyc(Rtj,Rtj_glm,params) ;
    
    % phase is computed modulo 2pi. in frequencies where the coherence is
    % not found to be statistically significant, the phase is set to NaN
    phi_sta_s = mod(phi_sta,2*pi) ;
    phi_sta_s(phi_sta_s>pi) = phi_sta_s(phi_sta_s>pi) - 2*pi ;
    phi_sta_s(Coh_sta<confC) = NaN ;

    phi_stc_s = mod(phi_stc,2*pi) ;
    phi_stc_s(phi_stc_s>pi) = phi_stc_s(phi_stc_s>pi) - 2*pi ;
    phi_stc_s(Coh_stc<confC) = NaN ;

    phi_glm_s = mod(phi_glm,2*pi) ;
    phi_glm_s(phi_glm_s>pi) = phi_glm_s(phi_glm_s>pi) - 2*pi ;
    phi_glm_s(Coh_glm<confC) = NaN ;

    save([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_validation_JK_' num2str(iJK) '.mat'],'confC','t','freq','SSt','SRt','SRt_sta','SRt_stc','SRt_glm','Coh_sta','Coh_stc','Coh_glm','phi_sta_s','phi_stc_s','phi_glm_s') ;
end
save([datadir 'Retina_cell_' num2str(icell) '_' stim_length{iL} '_loglikelihood.mat'],'LL') ;