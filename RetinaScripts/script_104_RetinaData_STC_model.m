% This script accompanies the Primer "Analysis of Neuronal Spike Trains, Deconstructed",
% by J. Aljadeff, B.J. Lansdell, A.L. Fairhall and D. Kleinfeld (2016) Neuron, 91 
% link to manuscript: http://dx.doi.org/10.1016/j.neuron.2016.05.039

% This script and script_104_RetinaData_STC_model.m compute the 
% Spike Triggered Covariance (STC) and the STC LN model (i.e. the 
% multi-dimensional linear nonlinear model where the dimesnions 
% or features - are the STA and STC) (Equations 18-24). The nonlinearity is 
% computed by averaging the response conditioned on the stimulus being in a
% specific bin (Equations 6 and 7). 

% The STC model is computed in five steps:
% (steps 1,2,3 in script_104_RetinaData_STC_significance.m and 
%  steps 4,5   in script_104_RetinaData_STC_model.m) 

% 1. computation of the spike triggered average. (Equations 14 and 15)
% 2. computation of the spike triggered covariance (Cs) and the matrix of 
%    covariance differences (dC) (Equations 18-22) 
% 3. finding significant STC dimensions by comparing to a null eigenvalue
%    distribution (Equation 23)
% 4. computation of the spiking nonlinearity (Equations 6 and 7)
% 5. building an LN for general stimulus by interpolating the spiking 
%    nonlinearity.

clear ;

workdir = uigetdir ; % select the parent directory through a GUI
cd(workdir) ;
datadir = 'RetinaData/' ;

if exist('histcn','file') == 0
    error('this script uses histcn.m, the multidimensional histogram function. please add it to the search path and rerun') ;
end

% maximal number of dimensions for model.
% if a number of STC features are found to be significant, it can be
% difficult to estimate the nonlinearity (Equations 6 and 7) because many
% stimulus bins may not be properly sampled. This sets a limit on the
% number of features that are used in the model. The code is written
% assuming kmax can be 1,2,3.
kmax = 2 ;         
                    
% a threshold for the projection of the STC features on the STA. 
% The STC features are computed independently of the STA and thus may have
% a large overlap with the STA. If this overlap is greater than the
% threshold chosen, we do not use the specific STC feature in the model.
% if the overlap is smaller than the threshold then the STC feature is
% used and the STA is projected out (Equation 24)
proj_th = 0.9 ;

% level of significance required to call an STC feature significant
% if an eigenvalue of dC is smaller than the a-th percentile of the null eigenvalue 
% distribution or larger than the 100-a percentile, it will be determined to be significant

a = 0.2 ;          
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
            load([datadir 'Retina_cell_' num2str(icell) '_sta_'       stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            load([datadir 'Retina_cell_' num2str(icell) '_stcsig_'    stim_length{iL} '_JK_' num2str(iJK) '.mat']) ;
            
            % 4. computation of the spiking nonlinearity 
            % ------------------------------------------
            
            % Finding STC features that have large overlap with STA and
            % removing them from significant STC features
            proj_stcsta = abs(vdC(:,isig)'*sta) ;
            isig(proj_stcsta > proj_th) = [] ; 
            
            % updating number of significant STC features
            nsig = length(isig) ;
        
            % the dimensionality of the model is set to 1 (STA) + number of
            % significant STC features (bounded by kmax)
            kmodel = 1 + min(nsig,kmax-1) ; 
            
            % stc is a matrix containing the significant STC features
            % here these features are ortho-normalized (Equation 24)
            stc = zeros(N,kmodel-1) ;
            for iisig = 1:kmodel-1
                stc(:,iisig) = vdC(:,isig(iisig)) - (vdC(:,isig(iisig))'*sta)*sta ;
                stc(:,iisig) = stc(:,iisig)/norm(stc(:,iisig)) ;
            end
            
            % f is the relevant stimulus subspace (STA + orthogonalized
            % significant STC modes)
            f = [sta stc] ;
            
            % projecting the stimulus onto the subspace f (Equation 5)
            % we assume that the response at time t depends on 
            % the stimulus only through z(t) (a K dimensional vector with 
            % K<<N).
            z = S*f ;
            
            % finding the range of z values to generate appropriate bins
            % (Equation 7)
            zx = max(abs(z(:)))+1e-3 ;
        
            bns = linspace(-zx,zx,nbns) ;  % bin edges
            dbns = bns(2)-bns(1) ;         % size of each bin (equal for all bins here for simplicity)
            ctr = bns(1:nbns-1)+dbns/2 ;   % bin centers
        
            bns_plot = linspace(-zx,zx,nbns_plot) ;
            dbns_plot = bns_plot(2)-bns_plot(1) ;
            ctr_plot = bns_plot(1:nbns_plot-1)+dbns_plot/2 ;
                                           % same for bins used for plotting 
            
            % average number of spikes per bin independent of stimulus
            nT = sum(R) ;
            
            
            switch kmodel
                % if K = 1 there are no significant STC features, so the
                % model is identical to the STA model
                case 1
                    % computing the nonlinear function g that converts stimulus 
                    % projections (z) to predicted firing rate (Equation 6). The 
                    % average response is computed at every bin of the projected 
                    % stimulus.
            
                    % Pf is the underlying stimulus distribution (projected on the 
                    % STA). normalized such that the integral over the distribution
                    % equals 1.
                    
                    g = zeros(nbns,1) ;
                    [Pf,~,~,bS] = histcn(z,bns) ;
                    for ib1 = 1:nbns
                        g(ib1) = mean(R(find(bS(:,1)==ib1))) ;
                    end
                    
                    Pf_plot  = histcn(z,bns_plot) ;    
                    Pf_plot  = Pf_plot/(T*dbns_plot) ; 
                
                % if K = 2 the model cell's response will depend on the projection 
                % of the stimulus on the STA and on a single STC feature
                case 2
                    
                    % computing the nonlinear function g that converts stimulus 
                    % projections (z) to predicted firing rate (Equation 6). The 
                    % average response is computed at every bin of the projected 
                    % stimulus. Also computing the marginal nonlinearities
                    % (Equation 10).
            
                    % Pf is the underlying stimulus distribution (projected on the 
                    % STA and STC). normalized such that the integral over the distribution
                    % equals 1. Also computed here are the marginals.
                    g = zeros(nbns,nbns) ;
                    g_a = zeros(nbns,1) ; % (STA marginal)
                    g_c = zeros(nbns,1) ; % (STC marginal)
                    
                    [Pf,~,~,bS] = histcn(z,bns,bns) ;
                    for ib1 = 1:nbns
                        g_a(ib1) = mean(R(bS(:,1)==ib1)) ;
                        g_c(ib1) = mean(R(bS(:,2)==ib1)) ;
                        for ib2 = 1:nbns
                            g(ib1,ib2) = mean(R(find((bS(:,1)==ib1).*(bS(:,2)==ib2)))) ;
                        end
                    end
                    

                    Pf_a  = histcn(z(:,1),bns) ;
                    Pf_a  = Pf_a/(T*dbns) ;             
                    
                    Pf_c  = histcn(z(:,2),bns) ;
                    Pf_c  = Pf_c/(T*dbns) ;
                    
                    Pf_a_plot  = histcn(z(:,1),bns_plot) ;
                    Pf_a_plot  = Pf_a_plot/(T*dbns_plot) ; 
                    
                    Pf_c_plot  = histcn(z(:,2),bns_plot) ;
                    Pf_c_plot  = Pf_c_plot/(T*dbns_plot) ;
                    
                % if K = 3 the model cell's response will depend on the projection 
                % of the stimulus on the STA and on a twp STC features
                case 3
                    % computing the nonlinear function g that converts stimulus 
                    % projections (z) to predicted firing rate (Equation 6). The 
                    % average response is computed at every bin of the projected 
                    % stimulus. Also computing the marginal nonlinearities
                    % (Equation 10).
                    
                    % Pf is the underlying stimulus distribution (projected on the 
                    % STA, STC_1, STC_2). normalized such that the integral over the distribution
                    % equals 1. Also computed here are the marginals.
                    g = zeros(nbns,nbns,nbns) ;
                    g_a  = zeros(nbns,1) ; % (STA   marginal)
                    g_c1 = zeros(nbns,1) ; % (STC_1 marginal)
                    g_c2 = zeros(nbns,1) ; % (STA_2 marginal)
                    [Pf,~,~,bS] = histcn(z,bns,bns,bns) ;
                    for ib1 = 1:nbns
                        g_a(ib1) = mean(R(bS(:,1)==ib1)) ;
                        g_c1(ib1) = mean(R(bS(:,2)==ib1)) ;
                        g_c2(ib1) = mean(R(bS(:,3)==ib1)) ;
                        for ib2 = 1:nbns
                            for ib3 = 1:nbns
                                g(ib1,ib2,ib3) = mean(R(find((bS(:,1)==ib1).*(bS(:,2)==ib2).*(bS(:,3)==ib3)))) ;
                            end
                        end
                    end
                    
                    Pf_a  = histcn(z(:,1),bns) ;
                    Pf_a  = Pf_a/(T*dbns) ; 
                    
                    Pf_c1  = histcn(z(:,2),bns) ;
                    Pf_c1  = Pf_c1/(T*dbns) ; 
                    
                    Pf_c2  = histcn(z(:,3),bns) ;
                    Pf_c2  = Pf_c2/(T*dbns) ; 
                    
                    Pf_a_plot  = histcn(z(:,1),bns_plot) ;
                    Pf_a_plot  = Pf_a_plot/(T*dbns_plot) ; 
                    
                    Pf_c1_plot  = histcn(z(:,2),bns_plot) ;
                    Pf_c1_plot  = Pf_c1_plot/(T*dbns_plot) ;
                    
                    Pf_c2_plot  = histcn(z(:,3),bns_plot) ;
                    Pf_c2_plot  = Pf_c2_plot/(T*dbns_plot) ;
            end
        
            g(isnan(g)) = 0 ;
            Pf = Pf/(T*dbns^kmodel) ;
            
            % 5. building an LN for general stimulus by interpolating the spiking nonlinearity
            % --------------------------------------------------------------------------------
            
            % the STC model is defined to be a 
            % linear interpolation of the spiking nonlinearity
            switch kmodel
                case 1 
                    stc_model = @(x)interpn(bns,g,x(:,1),'linear') ;
                    X = bns_plot ;
                case 2 
                    stc_model = @(x)interpn(bns,bns,g,x(:,1),x(:,2),'linear') ;
                    
                    [x1,x2] = ndgrid(bns_plot,bns_plot) ;
                    x1v = reshape(x1,nbns_plot^kmodel,1) ;
                    x2v = reshape(x2,nbns_plot^kmodel,1) ;
                    X = [x1v x2v] ;
                case 3
                    stc_model = @(x)interpn(bns,bns,bns,g,x(:,1),x(:,2),x(:,3),'linear') ;
                    [x1,x2,x3] = ndgrid(bns_plot,bns_plot,bns_plot) ;
                    x1v = reshape(x1,nbns_plot^kmodel,1) ;
                    x2v = reshape(x2,nbns_plot^kmodel,1) ;
                    x3v = reshape(x3,nbns_plot^kmodel,1) ;
                    X = [x1v x2v x3v] ;
            end
            
            % then the model is passed through a linear rectifier, to 
            % assure that the probability of a spike is non-negative
            % the infinitesimal constant offset (1e-8) ensures that the 
            % log-likelihood calculation does not diverge
            stc_model_rect = @(x) heaviside(stc_model(x)).*stc_model(x)+1e-8 ;
            stc_model_z = stc_model_rect(z) ;
            stc_model_z(isnan(stc_model_z)) = 0 ;

            % finally the model has to be renormalized such that the number 
            % of predicted spikes (on the training set) is fixed to what was 
            % found in the experiment. This step is needed because of the 
            % interpolation and rectification operations.  
            % The rectified interpolated model is scaled such that for the training 
            % set it predicts the total number of spikes found in the experiment.
            nT_model = sum(stc_model_z) ;
            stc_model_rect_norm = @(x) stc_model_rect(x)*nT/nT_model ;
            
            if kmodel > 1 
                stc_model_plot = reshape(stc_model_rect_norm(X),nbns_plot*ones(1,kmodel)) ;
            end
            switch kmodel
                case 1
                    save([datadir 'Retina_cell_' num2str(icell) '_stc_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'stc_model_rect_norm','f','kmodel','Pf','bns','ctr','bns_plot','ctr_plot') ;
                case 2
                    % for a two dimensional model we compute the Singular
                    % Value Decomposition of the STC model, to test whether
                    % or not the nonlinearity is separable (Equation 37)
                    stc_model_plot_svd = stc_model_plot ;
                    stc_model_plot_svd(isnan(stc_model_plot_svd)) = 0 ;
                    [stcnlin_U,stcnlin_S,stcnlin_V] = svd(stc_model_plot_svd) ;
                    save([datadir 'Retina_cell_' num2str(icell) '_stc_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'stc_model_rect_norm','f','kmodel','Pf','stc_model_plot','bns','ctr','bns_plot','ctr_plot','Pf_a','Pf_a_plot','Pf_c','Pf_c_plot','stcnlin_U','stcnlin_S','stcnlin_V','g_a','g_c') ;
                case 3
                    save([datadir 'Retina_cell_' num2str(icell) '_stc_' stim_length{iL} '_JK_' num2str(iJK) '.mat'],'stc_model_rect_norm','f','kmodel','Pf','stc_model_plot','bns','ctr','bns_plot','ctr_plot','Pf_a','Pf_a_plot','Pf_c1','Pf_c1_plot','Pf_c2','Pf_c2_plot','g_a','g_c1','g_c2') ; 
            end
        end
    end
end