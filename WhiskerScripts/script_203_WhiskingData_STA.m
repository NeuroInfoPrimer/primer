% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the Spike Triggered Average (STA) and the STA LN
%  model (i.e. the 1 dimensional linear nonlinear model where the single
%  dimension - or feature - is the STA). The nonlinearity is computed using
%  Bayes' rule as explained in the manuscript (Equations 8 and 9). 

% The STA model is computed in three steps: 
% 1. computation of the spike triggered avereage (i.e. reducing the 
%    dimensionality of the neuron's stimulus dependence to 1). (Equation
%    14)
% 2. computation of the spiking nonlinearity (Equations 8 and 9)
% 3. Building an LN model for a general stimulus by interpolating the 
%    spiking nonlinearity.

clear ; 
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

if exist('histcn','file') == 0
    error('this script uses histcn.m, the multidimensional histogram function. please add it to the search path and rerun') ;
end

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ; % Cell indices (used in file names)

sr = 1000 ;     % (Hz) sampling rate
N = 150 ;       % stimulus dimensionality 
ds = 2 ;        % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;    % delta t of stimulus 

nbns = 10 ;          % number of bins used to construct the stimulus distributions and estimate the nonlinearity
nbns_plot = 20 ;

tp = -dt*(1:1:N) ;   % the time vector associated with the stimulus history

nJK = 5 ;            % number of jack-knives 
for i = 1:7 
    % looping over all cells 
    
    for iJK = 1:nJK
    % loop over jack-knives. the STA model is estimated nJK = 5 times 
    
    load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ; 
    
    iR = R==1 ;        % binary response vector 
    
    S = S/N ;          % rescaled training stimulus. 
                       %  S is a (T x N) matrix meaning that 
                       %  the 1st dimension corresponds to time
                       %  the 2nd dimension corresponds to the stimulus 'pixels'
                       
    % 1. computation of the spike triggered avereage:
    % -----------------------------------------------
    
    sta = mean(S(iR,:),1) - mean(S,1) ; % implementation of Equations 14 and 15:
                                        % the 1st term on the RHS is the stimulus average conditioned on a spike
                                        % the 2nd term on the RHS is the unconditioned stimulus average 
                                        % here the STA is a (1 x N) vector
    
    sta = sta' / norm(sta) ;            % the STA is transposed to give a (N x 1) vector and normalized to 1 
    
    % 2. computation of the spiking nonlinearity 
    % ------------------------------------------
    
    z = S*sta ;                         % projecting the stimulus onto the STA (Equation 5)
                                        % we assume that the response at time t depends on 
                                        % the stimulus only through z(t) (a scalar quantity).
    
    zx = max(abs(z(:)))+1e-3 ;          % finding the range of z values to generate appropriate bins
                                        % (Equation 7)
        
    bns = linspace(-zx,zx,nbns) ;       % bin edges
    dbns = bns(2)-bns(1) ;              % size of each bin (equal for all bins here for simplicity)
    ctr = bns(1:nbns-1)+dbns/2 ;        % bin centers
        
    bns_plot = linspace(-zx,zx,nbns_plot) ;
    dbns_plot = bns_plot(2)-bns_plot(1) ;
    ctr_plot = bns_plot(1:nbns_plot-1)+dbns_plot/2 ;
                                        % same for bins used for plotting 

       
    Ps = mean(iR) ;                     % Probability of a spike independent of the stimulus
    
    Pf  = histcn(z,bns) ;               
    Pf  = Pf/(T*dbns) ;                 % Probability of stimulus (projected on STA feature) regardless of spiking
    
    Pfs = histcn(z(iR),bns) ;
    Pfs = Pfs/(T*Ps*dbns) ;             % Probability of stimulus (projected on STA feature) conditioned on a spike
    
    Psf = Pfs./Pf*Ps ;                  % Psf - probability of a spike given stimulus (projected on STA feature),
                                        % computed by invoking Bayes' rule (Equations 8 and 9). 
                                        % See explanations for the quantities Ps, Pf, Pfs in text following Equation 9.
    Psf(isnan(Psf)) = 0 ;               % setting NaNs to 0
        
    % same at higher resolution bins (for plotting)
    Pf_plot  = histcn(z,bns_plot) ;
    Pf_plot  = Pf_plot/(T*dbns_plot) ; 
    
    Pfs_plot = histcn(z(iR),bns_plot) ;
    Pfs_plot = Pfs_plot/(T*Ps*dbns_plot) ;
                
    Psf = Psf' ;                      % transposing Psf to give a (nbns x 1) vector
    
    % 3. Building an LN for general stimulus
    % --------------------------------------
        
    % first the STA model is defined to be a 
    % smooth interpolation of the spiking nonlinearity
    sta_model = @(x)interp1(ctr,Psf,x,'pchip') ; 
    
    % then the model is passed through a linear rectifier, to 
    % assure that the probability of a spike is non-negative
    % the infinitesimal constant offset (1e-8) ensures that the 
    % log-likelihood calculation does not diverge
    sta_model_rect = @(x) heaviside(sta_model(x)).*sta_model(x)+1e-8 ;  
    
    % finally the model has to be renormalized such that the number 
    % of predicted spikes (on the training set) is fixed to what was 
    % found in the experiment. This step is needed because of the 
    % smoothing and rectifying operations.  
    % The rectified smooth model is scaled such that for the training 
    % set it predicts the total number of spikes found in the experiment.
    Ps_model = sum(sta_model_rect(S*sta)) ;  
    sta_model_rect_norm = @(x) sta_model_rect(x)*Ps*T/Ps_model ;
    
    save([datadir 'VPM_cell_' Names{i} '_sta_JK_' num2str(iJK) '.mat'],'sta_model_rect_norm','sta','tp','N','ctr','ctr_plot','Pfs','Psf','Pf','Pf_plot','Pfs_plot') ; 
    
    end
end