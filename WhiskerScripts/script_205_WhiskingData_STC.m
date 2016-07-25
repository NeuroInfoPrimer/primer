% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the Spike Triggered Covariance (STC) and the STC LN
% model (i.e. the multi-dimensional linear nonlinear model where the dimesnions 
% - or features - are the STA and STC). The nonlinearity is computed using
% Bayes' rule as explained in the manuscript (Equations 8 and 9). 

% 1. computation of the spike triggered average. (Equations 14 and 15)
% 2. computation of the spike triggered covariance (Cs) and the matrix of 
%    covariance differences (dC) (Equations 18-22) 
% 3. finding significant STC dimensions by comparing to a null eigenvalue
%    distribution (Equation 23)
% 4. computation of the spiking nonlinearity (Equations 8 and 9)
% 5. building an LN for general stimulus by interpolating the spiking 
%    nonlinearity.

clear ; 
workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'WhiskerData/' ;

if exist('histcn','file') == 0
    error('this script uses histcn.m, the multidimensional histogram function. please add it to the search path and rerun') ;
end

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ; % Cell indices (used in file names)

sr = 1000 ;        % (Hz) sampling rate
N = 150 ;          % stimulus dimensionality 
ds = 2 ;           % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 

nJK = 5 ;          % number of jackknives

nspk = zeros(7,nJK) ; % a vector with the number of spikes for each cell 

tp = -dt*(1:1:N) ;
rep = 500 ;        % number of repeats used to compute the null distribution of covariance eigenvalues

a = 0.25 ;         % level of significance required to call an STC feature significant
                   % if an eigenvalue of dC is smaller than the a-th percentile of the null eigenvalue 
                   % distribution or larger than the 100-a percentile, it will be determined to be significant
nbns = 10 ;        % number of bins for each relevant stimulus dimension
                   % the stimulus frames (and correspoding spikes) will be 
                   % binned into nbins^2 bins 

nbns_plot = 25 ;

for i = 1:7
    % looping over all cells 

    for iJK = 1:nJK
    % loop over jack-knives.
        load([datadir 'VPM_cell_' Names{i} '_stim_resp_JK_' num2str(iJK) '.mat']) ; 
    
        nspk(i,iJK) = sum(R) ; % number of spikes 
        iR = R==1 ;            % binary response vector 
    
        S = S/N ;              % rescaled training stimulus. 
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

        % 2. computation of the spike triggered covariance:
        % -------------------------------------------------

        Cp = cov(S-repmat(mean(S,1),[T,1])) ;
                                            % implementation of Equations 18:
                                            % Cp is the covariance matrix of the mean subtracted underlying stimulus
 
        Cs = cov(S(iR,:)-repmat(mean(S(iR,:)),[nspk(i,iJK),1])) ; 
                                            % implementation of Equations 21:
                                            % Cs is the spike triggered covariance matrix
                                            
        % implementation of Equations 20, 22: 
        % the STC dimensions are found by diagonalizing dC - the matrix
        % of covariance differences
        dC = Cs-Cp ;               
        [vdC,edC] = eig(dC) ;
        [edC,idC] = sort(diag(edC)) ;
        vdC = vdC(:,idC) ;

        % 3. computation of the null eigenvalue distribution:
        % ---------------------------------------------------                                 

        % an array of covariance difference eigenvalues 
        % computed by assuming no association of the spike train
        % and the stimulus. (Equation 23)
        edCnull = zeros(1,N*rep) ;
        
        for ir = 1:rep                           
            % loop over number of repeats of shuffled spike trains 
        
            % random spike train obtained by shifting the real spike train by a random amount
            % with periodic boundary conditions. This preserves the possible burst structure and 
            % statistics of spike counts
            Rrnd = R(mod(randi(T)+(1:T),T)+1) ;  
            iRrnd = Rrnd == 1 ;

            % the null STC matrix is computed in the same way Cs is but conditioned on the 
            % random spike train rather than the real one (Equation 23)
            Crnd = cov(S(iRrnd,:)-repmat(mean(S(iRrnd,:)),[nspk(i,iJK),1])) ; 
            
            % N null eigenvalues are computed in each iteration, so in sum
            % the distributino of null eigenvalues will be constructed from 
            % rep x N eigenvalues.
            edCnull((ir-1)*N+(1:N)) = eig(Crnd-Cp) ;
        end
    
        max_enull = prctile(edCnull(:),100-a) ;  % upper bound of null distribution  
        min_enull = prctile(edCnull(:),a)  ;     % lower bound of null distribution
        vprct = linspace(a,100-a,N) ;
        prctenull = prctile(edCnull(:),vprct) ;  % percentiles of null distribution (all) 
    
        isig = [find(edC>max_enull) ; find(edC<min_enull)] ;
                                                 % the significant STC filters corrspond to eigenvalues of dC outside of the null distribution 
        nsig = length(isig) ;                    % number of significant STC features
        if nsig>0 
            [~,iis] = sort(abs(edC(isig)),'descend') ;
            isig = isig(iis) ;                   % orders the STC features according to the absolute value of their corresponding eigenvalues   
        end
    
        save([datadir 'VPM_cell_' Names{i} '_stcsig_JK_' num2str(iJK) '.mat'],'vdC','edC','edCnull','max_enull','min_enull','isig','nsig','a','vprct','prctenull') ;
    
        % If nsig > 0 we will construct an STC model. 
        % In the thalamic dataset the number of spikes/neuron ~ 1000, therefore the fit of the 
        % firing nonlinearity when the model dimension is > 2 is poor, so for the 
        % neurons with nsig > 0 we will construct 2D models (STA+STC1), only taking 
        % into account the most significant STC mode.
        
        if nsig > 0
            
            % stc is a matrix containing the significant STC features
            % here these features are ortho-normalized with respect to STA (Equation 24)
            stc = vdC(:,isig) - sta*(vdC(:,isig)'*sta)' ;
            for iisig = 1:nsig
                stc(:,iisig) = stc(:,iisig)/norm(stc(:,iisig)) ;
            end
            
            % the STC feature is the eigenvector of dC with largest eigenvalue (in absolute value)
        
            % projecting the stimulus onto the subspace spanned by STA and STC_1 (Equation 5)
            % we assume that the response at time t depends on 
            % the stimulus only through z(t) (a 2 dimensional vector)
            % z(t) is a [T x 2] array that holds the training stimulus 
            % projected on the STA (column 1) and on STC1 (column 2)
            z = zeros(T,2) ;                     
            z(:,1) = S*sta ;
            z(:,2) = S*stc(:,1) ;
        
            % finding the range of z values to generate appropriate bins
            % (Equation 7)
            zx = max(abs(z(:)))+1e-5 ;           
            
            % cell array of 2 bin vectors (here they are identical, but they do not have to be)
            bns = cell(1,2) ;                      
            bns{1} = linspace(-zx,zx,nbns) ;
            bns{2} = linspace(-zx,zx,nbns) ; 
        
            % the size of each bin (for computing proper probability distributions)   
            dbns = bns{1}(2)-bns{1}(1) ;         
            
            % same for high resolution bins used for plotting 
            bns_plot = cell(1,2) ;
            bns_plot{1} = linspace(-zx,zx,nbns_plot) ;
            bns_plot{2} = linspace(-zx,zx,nbns_plot) ; 
            dbns_plot = bns_plot{1}(2)-bns_plot{1}(1) ;
            
            % Ps  - probability of a spike
            %     - the number of spikes divided by the number of stimulus 'frames'
            Ps = nspk(i,iJK)/T ;                     
        
            % Pf  - probability of a stimulus projected on features
            %     - a normalized count of the projections of stimuli on the STA, STC1 in nbns^2 bins
            [Pf,ctr] = hist3(z,bns) ;Pf = Pf/(T*dbns^2) ;  
        
            % Pfs - probability of a stimulus projected on feature given a spike
            %     - a normalized count of the projections of stimuli on STA, STC1 in nbns^2 bins conditioned on a spike
            Pfs    = hist3(z(iR,:),bns) ;
            Pfs = Pfs/(nspk(i,iJK)*dbns^2) ;         
            
            % Psf - probability of a spike given stimulus projected on a feature
            %     - computed by invoking Bayes' rule (see details in text)
            Psf = Pfs./Pf*Ps ;
            Psf(isnan(Psf)) = 0 ;                % setting NaNs to 0
        
        
            % the same quantities, marginalized (Equation 10)
            % '_a' dependence on STA alone (averaged over STC dimension)
            % '_c' dependence on STC alone (averaged over STA dimension)
            Pf_a  = histcn(z(:,1),bns{1}) ;
            Pfs_a = histcn(z(iR,1),bns{1}) ;
            Pf_c  = histcn(z(:,2),bns{2}) ;
            Pfs_c = histcn(z(iR,2),bns{2}) ;
        
            Pf_a_plot  = histcn(z(:,1),bns_plot{1}) ;
            Pfs_a_plot = histcn(z(iR,1),bns_plot{1}) ;
            Pf_c_plot  = histcn(z(:,2),bns_plot{2}) ;
            Pfs_c_plot = histcn(z(iR,2),bns_plot{2}) ;
        
            Pf_a  = Pf_a/(T*dbns) ;
            Pfs_a = Pfs_a/(T*Ps*dbns) ;
            Pf_a_plot  = Pf_a_plot/(T*dbns_plot) ;  
            Pfs_a_plot = Pfs_a_plot/(T*Ps*dbns_plot) ;
            
            Pf_c  = Pf_c/(T*dbns) ;
            Pfs_c = Pfs_c/(T*Ps*dbns) ;
            Pf_c_plot  = Pf_c_plot/(T*dbns_plot) ;
            Pfs_c_plot = Pfs_c_plot/(T*Ps*dbns_plot) ;
            
            Psf_a = Pfs_a./Pf_a*Ps ;  
            Psf_a(isnan(Psf_a)) = 0 ;            %       setting NaNs to 0
            Psf_c = Pfs_c./Pf_c*Ps ;
            Psf_c(isnan(Psf_c)) = 0 ;            %       setting NaNs to 0

            Psf_a_plot = Pfs_a_plot./Pf_a_plot*Ps ;
            Psf_a_plot(isnan(Psf_a_plot)) = 0 ;  %       setting NaNs to 0
            Psf_c_plot = Pfs_c_plot./Pf_c_plot*Ps ;
            Psf_c_plot(isnan(Psf_c_plot)) = 0 ;  %       setting NaNs to 0

            % 5. building an LN for general stimulus by interpolating the spiking nonlinearity
            % --------------------------------------------------------------------------------

            % the STC model is defined to be a 
            % linear interpolation of the spiking nonlinearity
            stc_model = @(x)interpn(ctr{1},ctr{2},Psf,x(:,1),x(:,2),'linear') ;
    
            % then the model is passed through a linear rectifier, to 
            % assure that the probability of a spike is non-negative
            % the infinitesimal constant offset (1e-8) ensures that the 
            % log-likelihood calculation does not diverge
            stc_model_rect = @(x) heaviside(stc_model(x)).*stc_model(x)+1e-8 ;
        
            % finally the model has to be renormalized such that the number 
            % of predicted spikes (on the training set) is fixed to what was 
            % found in the experiment. This step is needed because of the 
            % interpolation and rectification operations.  
            % The rectified interpolated model is scaled such that for the training 
            % set it predicts the total number of spikes found in the experiment.
            Ps_model = sum(stc_model_rect(z)) ; 
            stc_model_rect_norm = @(x) stc_model_rect(x)*nspk(i,iJK)/Ps_model ;
        
            save([datadir 'VPM_cell_' Names{i} '_stc_JK_' num2str(iJK) '.mat'],'stc_model_rect_norm','sta','stc','tp','N','Pfs','Psf','Pf','nsig','isig','vdC','Pfs_a','Pfs_a_plot','Psf_a','Psf_a_plot','Pf_a','Pf_a_plot','Pfs_c','Pfs_c_plot','Psf_c','Psf_c_plot','Pf_c','Pf_c_plot','bns','nbns','dbns','bns_plot','nbns_plot','dbns_plot') ; 
        end
    end
end