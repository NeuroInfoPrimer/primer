% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the whisker phase tuning curve model (Equation 51)

clear ; 

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;

sr = 1000 ;               % (Hz) sampling rate
ds = 2 ;                  % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;              % delta t of stimulus 
nspk = zeros(1,7) ;       % a vector with the number of spikes for each cell 

nbns = 8 ;                % number of bins used to construct the stimulus distribution and estimate the nonlinearity

nJK = 5 ;                 % number of jack-knives 

for i = 1:7
    for iJK = 1:nJK 
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ; 
        nspk(i) = sum(R) ; % number of spikes 
        iR = R==1 ;        % binary response vector 
        
        % Ps  - probability of a spike                         
        %     - the number of spikes divided by the number of stimulus 'frames'    
        Ps = nspk(i)/T ;
        
        % Pp  - probability of a whisker phase
        %     - a normalized count of whisker phases in nbns bins
        [Pp,ctrs] = hist(Sph,nbns) ;
        
        dbns = ctrs(2)-ctrs(1) ;          % bin size
        Pp = Pp/T/dbns ;                  % normalizing Pp to 1     
    
        % Pps - probability of a whisker phase given a spike normalized to 1 
        Pps = hist(Sph(iR),ctrs)/nspk(i)/dbns ;
        
        % Psp - probability of a spike given whisker phase
        %     - computed by invoking Bayes' rule (Equation 51)
        Psp = Pps./Pp*Ps ;

        Psp(isnan(Psp)) = 0 ;             % setting NaNs to 0 
        Psp = Psp' ; 
    
        % first the tuning curve model is defined to be a
        % smooth interpolation of the binned onditioned probability distribution
        tun_model = @(x)interp1(ctrs,Psp,x,'pchip') ;
        
        % then the model passed through a linear rectifier, to
        % assure that the probability of a spike is non-negative
        % the infinitesimal constant offset (1e-8) ensures that the
        % loglikelihood calculation below does not diverge
        tun_model_rect = @(x) heaviside(tun_model(x)).*tun_model(x)+1e-8 ;
        
        % finally the model has to be renormalized such that the number of
        % of predicted spikes (on the training set) is fixed to what was
        % found in the experiment. This step is needed because of the
        % smoothing and rectifying operations
        Ps_model = sum(tun_model_rect(Sph)) ;
    
        % this normalization is done by computing the number of predicted
        % spikes from the smooth rectified model (Ps_model)                                                                 %
        % then the model multiplied by the ratio of the actual number of spikes and Ps_model                                                                        %
        tun_model_rect_norm = @(x) tun_model_rect(x)*nspk(i)/Ps_model ;
        
        save([datadir 'VPM_cell_' Names{i} '_tune_JK_' num2str(iJK) '.mat'],'ctrs','tun_model_rect_norm','Pp','Psp','Pps') ;
    end
end