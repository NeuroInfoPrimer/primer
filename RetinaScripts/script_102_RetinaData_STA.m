% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script computes the Spike Triggered Average (STA) and the STA LN
% model (i.e. the 1 dimensional linear nonlinear model where the single
% dimension - or feature - is the STA) (Equation 14). The nonlinearity is 
% computed by averaging the response conditioned on the stimulus being in a
% specific bin (Equations 6 and 7). 

% The STA model is computed in three steps: 
% 1. computation of the spike triggered avereage (i.e. reducing the 
%    dimensionality of the neuron's stimulus dependence to 1). (Equation
%    14)
% 2. computation of the spiking nonlinearity (Equations 6 and 7)
% 3. Building an LN model for a general stimulus by interpolating the 
%    spiking nonlinearity.

clear ; 

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

if exist('histcn','file') == 0
    error('this script uses histcn.m, the multidimensional histogram function. please add it to the search path and rerun') ;
end

stim_length = {'short2','short3','long'} ;

nbns = 15 ; % number of bins used to construct the underlying stimulus 
            % distribution and to estimate the nonlinearity 

nbns_plot = nbns*3 ; % number of bins used for plots

nJK = 5 ;            % number of jack-knives 

for icell = 3:3
        % results for cell 3 are presented in paper. to compute model for 
        % other cells change to a vector including all the cell
        % numbers you are interested in (from 1 to 53).
    for iL = 2:3
        % for the STA model we consider two stimulus configurations: 'long'
        % and 'short3'.

        for iJK = 1:nJK
            % loop over jack-knives. the STA model is estimated nJK = 5 times 
            load([datadir 'Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            
            % 1. computation of the spike triggered avereage:
            % -----------------------------------------------
    
            sta = S'*R/sum(R) - mean(S,1)' ;    % implementation of Equations 14 and 15:
                                                % the 1st term on the RHS is the stimulus average conditioned on a spike
                                                % the 2nd term on the RHS is the unconditioned stimulus average 
    
            sta = sta / norm(sta) ;             % the STA is normalized to 1 
        
    
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
        
            % computing the nonlinear function g that converts stimulus 
            % projections (z) to predicted firing rate (Equation 6). The 
            % average response is computed at every bin of the projected 
            % stimulus.
            
            % Pf is the underlying stimulus distribution (projected on the 
            % STA). normalized such that the integral over the distribution
            % equals 1.
            
            g = zeros(nbns,1) ;                 
            [Pf,~,~,bS] = histcn(z,bns) ;
            for ib = 1:nbns
                g(ib) = mean(R(bS==ib)) ;
            end
            g(isnan(g)) = 0 ;
            Pf  = Pf/(T*dbns) ;
            
            % same at higher resolution bins (for plotting)
            Pf_plot  = histcn(z,bns_plot) ;    
            Pf_plot  = Pf_plot/(T*dbns_plot) ; 
            
            % 3. Building an LN for general stimulus
            % --------------------------------------
            
            % first the STA model is defined to be a 
            % smooth interpolation of the spiking nonlinearity
            sta_model = @(x)interp1(bns,g,x,'pchip') ;                          
            
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
            nT_model = sum(sta_model_rect(z)) ; % the predicted total number of spikes before normalization
            sta_model_rect_norm = @(x) sta_model_rect(x)*sum(R)/nT_model ; 
            
            save([datadir 'Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'sta_model_rect_norm','sta','NT','NX','N','ctr','Pf','Pf_plot','ctr_plot','g','bns') ; 
        end
    end
end